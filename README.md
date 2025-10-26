# Discrete Elastic Rods

We implement a physically-based simulation of discrete elastic rods as described in [1]. This is done as part of the course "252-0546-00L Physically-Based Simulation in Computer Graphics" at ETH Zürich.

# Building

We use cmake's FetchContent functionality to retrieve our dependencies, however, some general dependencies should be installed system wide. This includes e.g. OpenGL.

There are two configure-presets available: One for debug builds and one for release builds.
For debug builds:
```sh
cmake --preset debug
cmake --build build_debug
./build_debug/discrete_elastic_rods
```
For release builds:
```sh
cmake --preset release
cmake --build build_release
./build_release/discrete_elastic_rods
```

A custom generator can of course be used as usual:
```sh
cmake --preset debug -GNinja
cmake --build build_debug
./build_debug/discrete_elastic_rods
```

# References

[1] Miklós Bergou et al. "Discrete elastic rods". In: ACM SIGGRAPH 2008 Papers. SIG-GRAPH '08. Los Angeles, California: Association for Computing Machinery, 2008. isbn: 9781450301121. doi: 10.1145/1399504.1360662. url: https://doi.org/10.1145/1399504.1360662.