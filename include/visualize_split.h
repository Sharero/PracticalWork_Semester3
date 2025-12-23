#ifndef VISUALIZE_SPLIT_H
#define VISUALIZE_SPLIT_H

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractEdges.h>
#include <vtkIntArray.h>
#include <vtkLookupTable.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkPyramid.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <map>
#include <vector>

class VTKVisualizer {
   public:
    static std::array<int, 4> orderPyramidBase(
        const std::array<int, 4>& base,
        const std::map<int, std::vector<double>>& coords);

    static std::array<int, 4> orderQuadPerPlane(
        const std::array<int, 4>& quad,
        const std::map<int, std::vector<double>>& coords);

    static void dumpUgridInfo(vtkUnstructuredGrid* ugrid);

    static vtkSmartPointer<vtkUnstructuredGrid> buildVTKGrid(
        const std::map<int, std::vector<double>>& global_node_coordinates,
        const std::vector<std::vector<int>>& pyramide_indicies,
        const std::vector<std::vector<int>>& tetrahedron_indicies);
};

#endif
