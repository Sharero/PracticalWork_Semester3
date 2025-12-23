#include "../include/visualize_split.h"

#include <algorithm>
#include <set>
#include <unordered_map>

std::array<int, 4> VTKVisualizer::orderPyramidBase(
    const std::array<int, 4>& base,
    const std::map<int, std::vector<double>>& coords) {
    double cx = 0;
    double cy = 0;

#pragma unroll 4
    for (int i : base) {
        cx += coords.at(i)[0];
        cy += coords.at(i)[1];
    }
    cx /= 4.0;
    cy /= 4.0;

    struct Item {
        int id;
        double ang;
    };

    std::array<Item, 4> items{};

#pragma unroll 4
    for (int i = 0; i < 4; ++i) {
        const auto& p = coords.at(base.at(i));
        items.at(i) = {.id = base.at(i),
                       .ang = std::atan2(p[1] - cy, p[0] - cx)};
    }

    std::ranges::sort(items, [](auto& a, auto& b) { return a.ang < b.ang; });

    return {items[0].id, items[1].id, items[2].id, items[3].id};
}

std::array<int, 4> VTKVisualizer::orderQuadPerPlane(
    const std::array<int, 4>& quad,
    const std::map<int, std::vector<double>>& coords) {
    std::array<double, 3> centroid{0, 0, 0};

#pragma unroll 4
    for (int i = 0; i < 4; ++i) {
        const auto& c = coords.at(quad.at(i));
        centroid[0] += c[0];
        centroid[1] += c[1];
        centroid[2] += c[2];
    }

#pragma unroll 4
    for (int k = 0; k < 3; ++k) {
        centroid.at(k) *= 0.25;
    }

    const auto& p0 = coords.at(quad[0]);
    const auto& p1 = coords.at(quad[1]);
    const auto& p2 = coords.at(quad[2]);

    std::array<double, 3> u{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    std::array<double, 3> v{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};

    std::array<double, 3> n{(u[1] * v[2]) - (u[2] * v[1]),
                            (u[2] * v[0]) - (u[0] * v[2]),
                            (u[0] * v[1]) - (u[1] * v[0])};

    double nn = std::sqrt((n[0] * n[0]) + (n[1] * n[1]) + (n[2] * n[2]));

    if (nn < 1e-12) {
        return orderPyramidBase(quad, coords);
    }

#pragma unroll 4
    for (int k = 0; k < 3; ++k) {
        n.at(k) /= nn;
    }

    std::array<double, 3> ux = u;
    double proj = (ux[0] * n[0]) + (ux[1] * n[1]) + (ux[2] * n[2]);
    ux = {ux[0] - (proj * n[0]), ux[1] - (proj * n[1]), ux[2] - (proj * n[2])};
    double ulen =
        std::sqrt((ux[0] * ux[0]) + (ux[1] * ux[1]) + (ux[2] * ux[2]));

    if (ulen < 1e-12) {
        return orderPyramidBase(quad, coords);
    }

#pragma unroll 4
    for (int k = 0; k < 3; ++k) {
        ux.at(k) /= ulen;
    }

    std::array<double, 3> uy = {(n[1] * ux[2]) - (n[2] * ux[1]),
                                (n[2] * ux[0]) - (n[0] * ux[2]),
                                (n[0] * ux[1]) - (n[1] * ux[0])};

    struct Item {
        int g;
        double ang;
        double r2;
    };

    std::array<Item, 4> items{};

#pragma unroll 4
    for (int i = 0; i < 4; ++i) {
        const auto& p = coords.at(quad.at(i));
        std::array<double, 3> r{p[0] - centroid[0], p[1] - centroid[1],
                                p[2] - centroid[2]};

        double uxv = (r[0] * ux[0]) + (r[1] * ux[1]) + (r[2] * ux[2]);
        double uyv = (r[0] * uy[0]) + (r[1] * uy[1]) + (r[2] * uy[2]);

        items.at(i) = {.g = quad.at(i),
                       .ang = std::atan2(uyv, uxv),
                       .r2 = (uxv * uxv) + (uyv * uyv)};
    }

    std::ranges::sort(items, [](const Item& a, const Item& b) {
        if (std::abs(a.ang - b.ang) > 1e-9) {
            return a.ang < b.ang;
        }
        return a.r2 < b.r2;
    });

    return {items[0].g, items[1].g, items[2].g, items[3].g};
}

void VTKVisualizer::dumpUgridInfo(vtkUnstructuredGrid* ugrid) {
    std::cout << "UGRID: points=" << ugrid->GetNumberOfPoints()
              << " cells=" << ugrid->GetNumberOfCells() << "\n";

    for (vtkIdType cid = 0; cid < ugrid->GetNumberOfCells(); ++cid) {
        vtkCell* cell = ugrid->GetCell(cid);
        int ctype = cell->GetCellType();
        std::cout << "Cell " << cid << " type=" << ctype
                  << " npts=" << cell->GetNumberOfPoints() << " pts: ";

#pragma unroll 4
        for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); ++i) {
            vtkIdType pid = cell->GetPointId(i);
            std::array<double, 3> p{};
            ugrid->GetPoint(pid, p.data());
            std::cout << pid << "(" << p[0] << "," << p[1] << "," << p[2]
                      << ") ";
        }
        std::cout << "\n";
    }
}

vtkSmartPointer<vtkUnstructuredGrid> VTKVisualizer::buildVTKGrid(
    const std::map<int, std::vector<double>>& global_node_coordinates,
    const std::vector<std::vector<int>>& pyramide_indicies,
    const std::vector<std::vector<int>>& tetrahedron_indicies) {
    std::set<int> used_indices;

    for (const auto& cell : pyramide_indicies) {
#pragma unroll 4
        for (int gi : cell) {
            used_indices.insert(gi);
        }
    }

    for (const auto& cell : tetrahedron_indicies) {
#pragma unroll 4
        for (int gi : cell) {
            used_indices.insert(gi);
        }
    }

    auto points = vtkSmartPointer<vtkPoints>::New();
    std::unordered_map<int, vtkIdType> global_to_vtk_id;
    global_to_vtk_id.reserve(used_indices.size());

#pragma unroll 4
    for (int gi : used_indices) {
        auto it = global_node_coordinates.find(gi);

        const std::vector<double>& c = it->second;
        vtkIdType pid = points->InsertNextPoint(c[0], c[1], c[2]);
        global_to_vtk_id[gi] = pid;
    }

    auto ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrid->SetPoints(points);

    auto elem_type = vtkSmartPointer<vtkIntArray>::New();
    elem_type->SetName("elemType");

    for (const auto& cell : pyramide_indicies) {
        auto pyr = vtkSmartPointer<vtkPyramid>::New();

        std::array<int, 4> base = {cell[0], cell[1], cell[2], cell[3]};

        auto ordered = orderQuadPerPlane(base, global_node_coordinates);

#pragma unroll 4
        for (int i = 0; i < 4; ++i) {
            pyr->GetPointIds()->SetId(i, global_to_vtk_id.at(ordered.at(i)));
        }

        pyr->GetPointIds()->SetId(4, global_to_vtk_id.at(cell[4]));

        ugrid->InsertNextCell(pyr->GetCellType(), pyr->GetPointIds());
        elem_type->InsertNextValue(1);
    }

    for (const auto& cell : tetrahedron_indicies) {
        auto tet = vtkSmartPointer<vtkTetra>::New();
#pragma unroll 4
        for (int i = 0; i < 4; ++i) {
            int g = cell[i];
            tet->GetPointIds()->SetId(i, global_to_vtk_id.at(g));
        }
        ugrid->InsertNextCell(tet->GetCellType(), tet->GetPointIds());
        elem_type->InsertNextValue(2);
    }

    ugrid->GetCellData()->AddArray(elem_type);
    ugrid->GetCellData()->SetActiveScalars("elemType");

    dumpUgridInfo(ugrid);
    return ugrid;
}