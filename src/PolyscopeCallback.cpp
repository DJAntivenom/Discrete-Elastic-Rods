/**
 * \file PolyscopeCallback.cpp
 * \brief Contains the definition of the polyscope callback.
 *
 * This file statically contains the variables used during simulation and
 * for driving the GUI.
 */

#include <DiscreteElasticRod.hpp>
#include <Eigen/Eigen>
#include <PolyscopeCallback.hpp>
#include <vector>

 /* Get rid of annoying unused ... warnings */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include <polyscope/curve_network.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#pragma GCC diagnostic pop

/**
 * \brief Container for all rods in the scene.
 *
 * Initially this is just one rod (unless we add the option
 * to start with multiple, like e.g. the chain-link in the
 * paper), and once a rod is cut, the number can grow.
 *
 * TODO: Think about if we really need this. We might also
 * get away with a single DiscreteElasticRod object which
 * can have multiple boundaries if we want to model a cut
 * rod. This might make collision detection simpler? But might
 * be less elegant.
 */
static std::vector<DiscreteElasticRod> rods;

/**
 * \brief Number of handles per rod.
 */
const int handles_count = 6;

/**
 * \brief The position of the handles to control the rods
 *
 * Block matrix of dimension (handles_count * rods.size()) x 3,
 * where each handles_count rows contains the handle positions
 * for a rod.
 */
static Eigen::MatrixXf handles;

/**
 * \brief The colors for the handles
 *
 * Must correspond to the dimension and content of handles.
 */
static Eigen::MatrixXf handle_colors;

/**
 * \brief Currently dragged handle. -1 for none.
 */
static int current_handle = -1;

/**
 * \brief Size of the timestep in seconds.
 */
static float delta_time = 0.001;

/**
 * @brief If true, the centerline will be drawn and the frame is transparent.
 */
static bool draw_centerline = true;

/**
 * @brief If true, the handles for controlling the rods will be drawn.
 */
static bool draw_handles = true;

/**
 * @brief Relative radius of the centerline.
 */
static float centerline_radius = 0.02;

/**
 * @brief Number of vertices to add.
 */
static int n_to_add = 10;

/**
 * @brief Maximal number of iterations newton method can run for in one timestep.
 */
static int max_newton_iterations = 5;

/**
 * \brief Written to by cout, see teebuffer below.
 */
static std::ostringstream log_stream;

/**
 * \brief Update the position of the handles for selected rod.
 */
static void updateRodHandles(int rod_index)
{
    const float direction_offset = 0.2f;
    const float rotation_offset = 0.1f;

    Eigen::MatrixX3f vertex_positions = rods[rod_index].getVertexPositions().reshaped(3, Eigen::AutoSize).transpose();

    Eigen::RowVector3f direction_handle_dir_start = -(vertex_positions.row(1) - vertex_positions.row(0)).normalized() * direction_offset;
    Eigen::RowVector3f orientation_handle_dir_start = Eigen::RowVector3f(-direction_handle_dir_start[1], direction_handle_dir_start[0], 0).normalized() * rotation_offset; // Any possible normal vector
    // TODO: Show a Kirchoff natural frame vector instead for the orientation

    uint64_t vertex_count = vertex_positions.rows();
    Eigen::RowVector3f direction_handle_dir_end = -(vertex_positions.row(vertex_count - 1) - vertex_positions.row(vertex_count - 2)).normalized() * direction_offset;
    Eigen::RowVector3f orientation_handle_dir_end = Eigen::RowVector3f(-direction_handle_dir_end[1], direction_handle_dir_end[0], 0).normalized() * rotation_offset;

    int i = 0;
    // Start position handle, head of the rod
    handles(rod_index + i++, Eigen::seq(0, 2)) = vertex_positions.row(0);
    // Start direction handle, offset from the head
    handles(rod_index + i++, Eigen::seq(0, 2)) = vertex_positions.row(0) + direction_handle_dir_start;
    // Rotation setting, offset to the side
    handles(rod_index + i++, Eigen::seq(0, 2)) = vertex_positions.row(0) + orientation_handle_dir_start;

    int verts = vertex_positions.rows();
    // End position handle, tail of the rod
    handles(rod_index + i++, Eigen::seq(0, 2)) = vertex_positions.row(verts - 1);
    // End direction handle, offset from the tail
    handles(rod_index + i++, Eigen::seq(0, 2)) = vertex_positions.row(verts - 1) - direction_handle_dir_end;
    // Rotation setting, offset to the side
    handles(rod_index + i++, Eigen::seq(0, 2)) = vertex_positions.row(verts - 1) + orientation_handle_dir_end;
}

/**
 * \brief Initialize a new rod.
 *
 * TODO: Where to put these helper functions?
 */
static void initializeRod(int n_vertices)
{
    rods.emplace_back(n_vertices);

    // Add the new handle positions
    handles.resize(handles.rows() + handles_count, 3);
    handle_colors.resize(handle_colors.rows() + handles_count, 3);

    // Initialize positions
    const uint64_t rod_count = rods.size();
    for (uint64_t rod_index = 0; rod_index < rod_count; ++rod_index)
    {
        updateRodHandles(rod_index);
    }

    // Set the colors
    handle_colors(Eigen::seq(handles.rows() - handles_count, handles.rows() - 1), Eigen::seq(0, 2)).setZero();
    // Blue for the endpoints
    handle_colors(handles.rows() - 6, Eigen::seq(0, 2)) = handle_colors(handles.rows() - 3, Eigen::seq(0, 2)) = Eigen::Vector3f(0.0, 0.0, 1.0).transpose();
    // Red for the direction
    handle_colors(handles.rows() - 5, Eigen::seq(0, 2)) = handle_colors(handles.rows() - 2, Eigen::seq(0, 2)) = Eigen::Vector3f(1.0, 0.0, 0.0).transpose();
    // Green for the orientation
    handle_colors(handles.rows() - 4, Eigen::seq(0, 2)) = handle_colors(handles.rows() - 1, Eigen::seq(0, 2)) = Eigen::Vector3f(0.0, 1.0, 0.0).transpose();
}

/**
 * \brief Tries applying handle update to the rod
 *
 * TODO: Should this apply a force?
 */
static void applyHandleUpdate(uint64_t handle_index, const Eigen::Vector3f &new_position)
{
    uint64_t handle_type = handle_index % handles_count;
    uint64_t rod_index = handle_index / handles_count;

    // TODO: Add other handles
    switch (handle_type)
    {
    case 0:
        rods[rod_index].setVertexPosition(0, new_position);
        updateRodHandles(rod_index);
        break;
    case 3:
        rods[rod_index].setVertexPosition(rods[rod_index].getVertexPositions().cols() - 1, new_position);
        updateRodHandles(rod_index);
        break;
    default:
        break;
    }
}

static Eigen::Vector3f findNearestPointToRay(const Eigen::Vector3f &origin, const Eigen::Vector3f &direction, const Eigen::Vector3f &point)
{
    float t = (point - origin).dot(direction);
    return origin + t * direction;
}

/**
 * \brief Convert screen coordinates to world position
 *
 * Source: polyscope::view::screenCoordsToWorldRay
 * Shortened to return its intermediate value, `worldPos`.
 */
Eigen::Vector3f screenCoordsToWorldPos(glm::vec2 screenCoords) {
    auto view = polyscope::view::getCameraViewMatrix();
    auto proj = polyscope::view::getCameraPerspectiveMatrix();
    glm::vec4 viewport = { 0., 0., polyscope::view::windowWidth, polyscope::view::windowHeight };

    glm::vec3 screenPos3{ screenCoords.x, polyscope::view::windowHeight - screenCoords.y, 0. };
    auto worldPos = glm::unProject(screenPos3, view, proj, viewport);
    return Eigen::Vector3f(worldPos.x, worldPos.y, worldPos.z);
}

/**
 * \brief Creates the GUI holding the controls of the simulation.
 */
static void makeConfigWindow()
{
    if (ImGui::CollapsingHeader("Simulation settings", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("Time-step", &delta_time, 0.001);
        ImGui::InputInt("Max number of newton iterations", &max_newton_iterations);

        if (rods.size() == 0)
        {
            ImGui::InputInt("Vertices n", &n_to_add);
            if (ImGui::Button("Create Rod"))
            {
                initializeRod(n_to_add);
            }
        }
    }

    if (ImGui::CollapsingHeader("Visualization settings", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Draw handles", &draw_handles);
        ImGui::Checkbox("Draw centerline", &draw_centerline);
        if (draw_centerline)
        {
            ImGui::DragFloat("Relative Radius", &centerline_radius, 0.01f, 0.f, 1.f);
        }

        ImGui::Checkbox("Automatic bounding-box", &polyscope::options::automaticallyComputeSceneExtents);
    }

    if (ImGui::CollapsingHeader("Debug settings", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Print debug stuff"))
        {
            if (rods.size() > 0)
            {
                std::cout << rods[0].getEdges() << "\n";
            }
        }
    }
}

/**
 * \brief Creates the GUI holding information for the user.
 */
static void makeAnalysisWindow()
{
    auto lower_corner = std::get<0>(polyscope::state::boundingBox);
    auto upper_corner = std::get<1>(polyscope::state::boundingBox);
    ImGui::LabelText("Scene extent", "(%.1f,%.1f,%.1f)x(%.1f,%.1f,%.1f)",
                     lower_corner.x, lower_corner.y, lower_corner.z,
                     upper_corner.x, upper_corner.y, upper_corner.z);

    ImGui::TextWrapped("Log:\n%s", log_stream.str().c_str());
}

static void updateViewerData()
{
    /* If the user requests it, calculate scene extents */
    if (!polyscope::options::automaticallyComputeSceneExtents)
    {
        polyscope::state::boundingBox = std::make_tuple(glm::vec3(-1, -1, -1),
                                                        glm::vec3(1, 1, 1));
    }

    polyscope::removeAllStructures();

    if (draw_centerline)
    {
        const uint64_t rod_count = rods.size();
        for (uint64_t rod_index = 0; rod_index < rod_count; ++rod_index)
        {
            Eigen::MatrixX3f vertex_positions = rods[rod_index].getVertexPositions().transpose();

            auto lines = polyscope::registerCurveNetworkLine("Centerline_" + std::to_string(rod_index), vertex_positions);

            lines->setRadius(centerline_radius, false);
        }
    }

    const double handle_radius = 0.03;
    const std::string handle_structure_name = "Handles";

    if (draw_handles && !rods.empty())
    {
        auto pointcloud = polyscope::registerPointCloud(handle_structure_name, handles);
        pointcloud->setPointRadius(handle_radius, false);
        auto pointcloud_color = pointcloud->addColorQuantity("Color", handle_colors);
        pointcloud_color->setEnabled(true);
    }

    if (!ImGui::GetIO().WantCaptureMouse)
    {
        // Move currently selected handle accoring to the mouse delta
        if (current_handle != -1)
        {
            // Camera and Mouse parameters
            ImVec2 imgui_mouse_screen_pos = ImGui::GetIO().MousePos;
            glm::vec2 poly_mouse_screen_pos(imgui_mouse_screen_pos.x, imgui_mouse_screen_pos.y);
            ImVec2 imgui_mouse_screen_pos_prev = ImGui::GetIO().MousePosPrev;
            glm::vec2 poly_mouse_screen_pos_prev(imgui_mouse_screen_pos_prev.x, imgui_mouse_screen_pos_prev.y);

            Eigen::Vector3f mouse_position = screenCoordsToWorldPos(poly_mouse_screen_pos);
            Eigen::Vector3f mouse_position_prev = screenCoordsToWorldPos(poly_mouse_screen_pos_prev);

            glm::vec3 poly_mouse_direction = polyscope::view::screenCoordsToWorldRay(poly_mouse_screen_pos);
            Eigen::Vector3f mouse_direction(poly_mouse_direction.x, poly_mouse_direction.y, poly_mouse_direction.z);
            glm::vec3 poly_mouse_direction_prev = polyscope::view::screenCoordsToWorldRay(poly_mouse_screen_pos_prev);
            Eigen::Vector3f mouse_direction_prev(poly_mouse_direction_prev.x, poly_mouse_direction_prev.y, poly_mouse_direction_prev.z);

            auto handle = handles.row(current_handle);

            // Find the distance by how much to move the point
            // First, find the nearest point to the handle on the current and previous mouse ray. Then, move the handle by that distance.

            Eigen::Vector3f mouse_closest_point = findNearestPointToRay(mouse_position, mouse_direction, handle);
            Eigen::Vector3f mouse_closest_point_prev = findNearestPointToRay(mouse_position_prev, mouse_direction_prev, handle);

            // Raw mouse delta
            Eigen::Vector3f mouse_delta = mouse_closest_point - mouse_closest_point_prev;
            // Normalize this to be perpendicular to the current mouse direction
            // TODO: Validate that this feels good to move stuff with
            mouse_delta = mouse_delta - mouse_delta.dot(mouse_direction) * mouse_direction;

            Eigen::Vector3f new_position = handle + (mouse_delta).transpose();

            applyHandleUpdate(current_handle, new_position);
        }

        // Check if we are clicking to select a handle
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left, false))
        {
            // Clear out old handle selection, if exists
            current_handle = -1;

            // Find nearest handle
            ImVec2 imgui_mouse_pos = ImGui::GetIO().MousePos;
            polyscope::PickResult pick_result = polyscope::pickAtScreenCoords(glm::vec2(imgui_mouse_pos.x, imgui_mouse_pos.y));
            if (pick_result.isHit && pick_result.structureName == handle_structure_name)
            {
                Eigen::Vector3f hit_location(pick_result.position.x, pick_result.position.y, pick_result.position.z);

                // Look for handle at the same position, since we cannot tell from the pick result directly.
                const uint64_t handles_size = handles.rows();
                for (uint64_t handle_index = 0; handle_index < handles_size; handle_index++)
                {
                    float distance = (hit_location - handles.row(handle_index).transpose()).norm();
                    // Note: Distance is typically 1.3ish times the handle radius. 2x should detect but not generate false-positives
                    // No idea why it's larger than the radius.
                    if (distance < 2 * handle_radius)
                    {
                        current_handle = handle_index;
                        // Disable camera rotation
                        polyscope::state::doDefaultMouseInteraction = false;
                        break;
                    }
                }
            }
        }

        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left))
        {
            // Set handle not controlled
            current_handle = -1;
            // Re-enable camera rotation
            polyscope::state::doDefaultMouseInteraction = true;
        }
    }
}

void polyscopeCallback()
{
    /// Main menu window on screen left.
    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Once);
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x * 0.25, ImGui::GetIO().DisplaySize.y),
                             ImGuiCond_Once);
    ImGui::Begin("Menu");
    makeConfigWindow();
    ImGui::End();

    /// Additional "Analysis" window on right, content provided by SubApp.
    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x * 0.75, 0), ImGuiCond_Once);
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x * 0.25, ImGui::GetIO().DisplaySize.y),
                             ImGuiCond_Once);
    ImGui::Begin("Analysis");
    makeAnalysisWindow();
    ImGui::End();

    /// Update display.
    updateViewerData();

    /// Main update of rods
    for (auto &rod : rods)
    {
        rod.update(delta_time, static_cast<size_t>(max_newton_iterations));
    }
}

/////////////////////////////////////////////////////
////////// COUT CAPTURING IMPLEMENTATION ////////////
// https://wordaligned.org/articles/cpp-streambufs //
/////////////////////////////////////////////////////

#include <streambuf>
class teebuf : public std::streambuf
{
public:
    // Construct a streambuf which tees output to both input
    // streambufs.
    teebuf(std::streambuf *sb1) : sb1(sb1), sb2(log_stream.rdbuf()) {}
    teebuf() = default;

private:
    // This tee buffer has no buffer. So every character "overflows"
    // and can be put directly into the teed buffers.
    virtual int overflow(int c)
    {
        if (c == EOF)
        {
            return !EOF;
        }
        else
        {
            int const r1 = sb1->sputc(c);
            int const r2 = sb2->sputc(c);
            return r1 == EOF || r2 == EOF ? EOF : c;
        }
    }

    // Sync both teed buffers.
    virtual int sync()
    {
        int const r1 = sb1->pubsync();
        int const r2 = sb2->pubsync();
        return r1 == 0 && r2 == 0 ? 0 : -1;
    }

private:
    std::streambuf *sb1;
    std::streambuf *sb2;
};

static teebuf gui_tee_buffer;

void overrideCoutBuffer()
{
    gui_tee_buffer = teebuf(std::cout.rdbuf());
    std::cout.rdbuf(&gui_tee_buffer);
}
