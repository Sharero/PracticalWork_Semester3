#ifndef GRID_H
#define GRID_H

#include <map>
#include <vector>

#include "../include/visualize_split.h"

class Grid {
   private:
    VTKVisualizer visualizer;

    std::vector<std::vector<int>> tetrahedron_indicies;

    std::vector<std::vector<int>> pyramide_indicies;

   public:
    static std::vector<int> buildLocalMappingForFixedFace(
        const std::vector<int>& global_nodes_indicies,
        const std::map<int, std::vector<double>>& global_nodes_coordinates,
        const std::vector<int>& fixed_base_global_indicies);

    static std::vector<std::vector<int>> applyTemplate(
        const std::vector<int>& local_mapping);

    static void ensurePositiveTets(
        std::vector<std::vector<int>>& elems,
        const std::map<int, std::vector<double>>& global_nodes_coordinates);

    void splitHexAnyFixedFace(
        const std::vector<int>& global_nodes_indicies,
        const std::map<int, std::vector<double>>& global_nodes_coordinates,
        const std::vector<int>& fixed_base_global_indicies);

    void printElements();

    void visualizeSplittedGrid(
        const std::map<int, std::vector<double>>& global_nodes_coordinates);
};

#endif