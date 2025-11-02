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

 /* Get rid of anyoing unused ... warinings */
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
const int handles_count = 4;

/**
 * \brief The position of the handles to control the rods
 *
 * Block matrix of dimension (handles_count * rods.size()) x 3,
 * where each handles_count rows contains the handle positions
 * for a rod.
 */
static Eigen::MatrixXf handles;

/**
 * \brief The colors for the handles. Must correspond to the
 * dimension of handles.
 */
static Eigen::MatrixXf handle_colors;

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
 * \brief Written to by cout, see teebuffer below.
 */
static std::ostringstream log_stream;

/**
 * \brief Update the position of the handles for selected rod.
 */
static void update_rod_handles(int rod_index)
{

    // TODO: Make point sizes, distances, colors into a const?
    // TODO: This is currently inverted, so that the handles point in and don't destroy the scaling. Fix this after asking a TA.
    const float direction_offset = -0.2f;

    Eigen::MatrixX3f vertex_positions = rods[rod_index].getVertexPositions().reshaped(3, Eigen::AutoSize).transpose();

    // TODO: Cleanup duplication here, initialize rotation handles too
    // Handle 1: Start position handle, head of the rod
    handles(rod_index, Eigen::seq(0, 2)) = vertex_positions.row(0);
    // Handle 2: Start direction handle, offset from the head
    handles(rod_index + 1, Eigen::seq(0, 2)) = vertex_positions.row(0) - (vertex_positions.row(1) - vertex_positions.row(0)).normalized() * direction_offset;

    int verts = vertex_positions.rows();
    // Handle 3: End position handle, tail of the rod
    handles(rod_index + 2, Eigen::seq(0, 2)) = vertex_positions.row(verts - 1);
    // Handle 4: End direction handle, offset from the tail
    handles(rod_index + 3, Eigen::seq(0, 2)) = vertex_positions.row(verts - 1) - (vertex_positions.row(verts - 2) - vertex_positions.row(verts - 1)).normalized() * direction_offset;
}

/**
 * \brief Initialize a new rod.
 *
 * TODO: Where to put these helper functions?
 */
static void initialize_rod(int n_vertices)
{
    rods.emplace_back(n_vertices);

    // Add the new handle positions
    handles.resize(handles.rows() + handles_count, 3);
    handle_colors.resize(handle_colors.rows() + handles_count, 3);

    // Initialize positions
    const uint64_t rod_count = rods.size();
    for (uint64_t rod_index = 0; rod_index < rod_count; ++rod_index)
    {
        update_rod_handles(rod_index);
    }

    // Set the colors
    handle_colors(Eigen::seq(handles.rows() - handles_count, handles.rows() - 1), Eigen::seq(0, 2)).setZero();
    // Blue for the endpoints
    handle_colors(handles.rows() - 4, Eigen::seq(0, 2)) = handle_colors(handles.rows() - 2, Eigen::seq(0, 2)) = Eigen::Vector3f(0.0, 0.0, 1.0).transpose();
    // Red for the direction
    handle_colors(handles.rows() - 3, Eigen::seq(0, 2)) = handle_colors(handles.rows() - 1, Eigen::seq(0, 2)) = Eigen::Vector3f(1.0, 0.0, 0.0).transpose();
}

/**
 * \brief Creates the GUI holding the controls of the simulation.
 */
static void makeConfigWindow()
{
    ImGui::InputFloat("Time-step", &delta_time, 0.001);
    ImGui::Checkbox("Draw handles", &draw_handles);
    ImGui::Checkbox("Draw centerline", &draw_centerline);
    if (draw_centerline)
    {
        ImGui::DragFloat("Relative Radius", &centerline_radius, 0.01f, 0.f, 1.f);
    }
    if (rods.size() == 0)
    {
        ImGui::InputInt("Vertices n", &n_to_add);
        if (ImGui::Button("Create Rod"))
        {
            initialize_rod(n_to_add);
        }
    }
    ImGui::Checkbox("Automatic bounding-box", &polyscope::options::automaticallyComputeSceneExtents);
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
    // auto points = polyscope::registerPointCloud("Points" + std::to_string(i), curr_data.points);
    // auto color = points->addColorQuantity("Color", curr_data.points_c);
    // color->setEnabled(true);

    // points->setPointRadius(curr_data.point_size);
    // if (curr_data.points_quad)
    // {
    //     points->setPointRenderMode(polyscope::PointRenderMode::Quad);
    // }
    // points->setTransparency(curr_data.alpha);

    if (draw_centerline)
    {
        const uint64_t rod_count = rods.size();
        for (uint64_t rod_index = 0; rod_index < rod_count; ++rod_index)
        {
            Eigen::MatrixX3f vertex_positions = rods[rod_index].getVertexPositions().reshaped(3, Eigen::AutoSize).transpose();

            auto lines = polyscope::registerCurveNetworkLine("Centerline_" + std::to_string(rod_index), vertex_positions);

            lines->setRadius(centerline_radius);
        }
    }

    if (draw_handles && !rods.empty())
    {
        auto pointcloud = polyscope::registerPointCloud("Handles", handles);
        pointcloud->setPointRadius(0.03); // TODO: Make this a constant
        auto pointcloud_color = pointcloud->addColorQuantity("Color", handle_colors);
        pointcloud_color->setEnabled(true);

        // std::cout << "Points: " << points << std::endl;
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
        rod.update(delta_time);
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
