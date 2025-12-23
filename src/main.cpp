#include <array>
#include <fstream>
#include <iostream>
#include <vector>

#include "../include/grid.h"

constexpr auto GLOBAL_NODE_INDICIES_INPUT_FILE_NAME =
    std::to_array("./tests/5_test/global_node_indicies.txt");

constexpr auto GLOBAL_NODE_COORDINATES_INPUT_FILE_NAME =
    std::to_array("./tests/5_test/global_node_coordinates.txt");

int main() {
    Grid grid;

    std::ifstream global_indicies_input_file(
        GLOBAL_NODE_INDICIES_INPUT_FILE_NAME.data());

    int count_node_indicies = 0;
    int count_node_coordinates = 0;

    global_indicies_input_file >> count_node_indicies;
    std::vector<std::array<int, 8>> global_nodes_indicies(count_node_indicies /
                                                          8);

    for (int i = 0; i < count_node_indicies / 8; i++) {
        std::array<int, 8> temp{};
#pragma unroll 4
        for (int j = 0; j < 8; j++) {
            global_indicies_input_file >> temp.at(j);
        }
        global_nodes_indicies[i] = temp;
    }

    global_indicies_input_file.close();

    std::ifstream global_coordinates_input_file(
        GLOBAL_NODE_COORDINATES_INPUT_FILE_NAME.data());

    global_coordinates_input_file >> count_node_coordinates;

    std::map<int, std::vector<double>> global_node_coordinates;

#pragma unroll 4
    for (int i = 0; i < count_node_coordinates; i++) {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        global_coordinates_input_file >> x >> y >> z;
        global_node_coordinates[i] = {x, y, z};
    }

    global_coordinates_input_file.close();

    std::vector<int> fixed_face_local_indicies = {0, 1, 2, 3};

#pragma unroll 4
    for (auto& global_nodes_index : global_nodes_indicies) {
        std::vector<int> nodes8(8);

#pragma unroll 4
        for (int k = 0; k < 8; ++k) {
            nodes8[k] = global_nodes_index.at(k);
        }

        std::vector<int> fixed_face_local_indicies = {nodes8[0], nodes8[1],
                                                      nodes8[2], nodes8[3]};

        grid.splitHexAnyFixedFace(nodes8, global_node_coordinates,
                                  fixed_face_local_indicies);
    }

    // grid.printElements();
    grid.visualizeSplittedGrid(global_node_coordinates);
    return 0;
}
