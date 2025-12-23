#include "../include/grid.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>

#include "../include/math_utils.h"

static constexpr std::array<int, 5> PYRAMID_TEMPLATE = {0, 1, 2, 3, 4};

static constexpr std::array<std::array<int, 4>, 4> TETRAHEDRON_TEMPLATE{
    {{{1, 3, 4, 5}}, {{2, 3, 4, 6}}, {{3, 4, 5, 6}}, {{3, 5, 6, 7}}}};

constexpr double EPS = 1e-9;

std::vector<int> Grid::buildLocalMappingForFixedFace(
    const std::vector<int>& global_nodes_indicies,
    const std::map<int, std::vector<double>>& global_nodes_coordinates,
    const std::vector<int>& fixed_base_global_indicies) {
    std::vector<std::pair<int, std::vector<double>>> base;
    base.reserve(4);

    std::set<int> base_set(fixed_base_global_indicies.begin(),
                           fixed_base_global_indicies.end());

#pragma unroll 4
    for (int g : fixed_base_global_indicies) {
        auto it = global_nodes_coordinates.find(g);
        base.emplace_back(g, it->second);
    }

    std::vector<double> a = base[0].second;
    std::vector<double> b = base[1].second;
    std::vector<double> c = base[2].second;

    std::vector<double> n = MathHelper::calculateCrossProductOfVectors(
        MathHelper::subtractVectors(b, a), MathHelper::subtractVectors(c, a));

    if (MathHelper::calculateNormOfVector(n) < EPS) {
        c = base[3].second;
        n = MathHelper::calculateCrossProductOfVectors(
            MathHelper::subtractVectors(b, a),
            MathHelper::subtractVectors(c, a));
        if (MathHelper::calculateNormOfVector(n) < EPS) {
            n = {0, 0, 0};
#pragma unroll 4
            for (int i = 0; i < 4; ++i) {
                std::vector<double> p0 = base[i].second;
                std::vector<double> p1 = base[(i + 1) % 4].second;
                std::vector<double> p2 = base[(i + 2) % 4].second;
                std::vector<double> c =
                    MathHelper::calculateCrossProductOfVectors(
                        MathHelper::subtractVectors(p1, p0),
                        MathHelper::subtractVectors(p2, p0));
                n = MathHelper::additionVectors(n, c);
            }
        }
    }
    n = MathHelper::calculateNormalizedVector(n);

    std::vector<double> u =
        MathHelper::subtractVectors(base[1].second, base[0].second);
    if (MathHelper::calculateNormOfVector(u) < EPS) {
        u = MathHelper::subtractVectors(base[2].second, base[0].second);
    }
    u = MathHelper::calculateNormalizedVector(u);
    double proj = MathHelper::calculateScalarProductOfVectors(u, n);
    u = MathHelper::subtractVectors(
        u, MathHelper::multiplyVectorByMultiplier(n, proj));
    u = MathHelper::calculateNormalizedVector(u);
    std::vector<double> v = MathHelper::calculateCrossProductOfVectors(n, u);
    v = MathHelper::calculateNormalizedVector(v);

    std::vector<double> centroid = {0, 0, 0};
#pragma unroll 4
    for (auto& p : base) {
        centroid = MathHelper::additionVectors(centroid, p.second);
    }
    centroid = MathHelper::multiplyVectorByMultiplier(centroid, 1.0 / 4.0);

    struct Item {
        int g;
        double uc;
        double vc;
        std::vector<double> pos;
    };
    std::vector<Item> items;
    items.reserve(4);
#pragma unroll 4
    for (auto& p : base) {
        std::vector<double> r = MathHelper::subtractVectors(p.second, centroid);
        double uc = MathHelper::calculateScalarProductOfVectors(r, u);
        double vc = MathHelper::calculateScalarProductOfVectors(r, v);
        items.push_back({p.first, uc, vc, p.second});
    }

    std::ranges::sort(items, [](const Item& a, const Item& b) {
        if (std::abs(a.vc - b.vc) > 1e-8) {
            return a.vc < b.vc;
        }
        return a.uc < b.uc;
    });

    std::vector<int> others;
    others.reserve(4);
#pragma unroll 4
    for (int i = 0; i < 8; ++i) {
        int g = global_nodes_indicies[i];
        if (!base_set.contains(g)) {
            others.push_back(g);
        }
    }

    auto find_opposite = [&](const std::vector<double>& base_pos) -> int {
        double best = -1e308;
        int bestg = -1;
#pragma unroll 4
        for (int g : others) {
            std::vector<double> r = MathHelper::subtractVectors(
                global_nodes_coordinates.at(g), base_pos);
            double val = MathHelper::calculateScalarProductOfVectors(r, n);
            if (val > best) {
                best = val;
                bestg = g;
            }
        }
        return bestg;
    };

    std::vector<int> used_opp(others.size(), 0);
    std::vector<int> local(8);
#pragma unroll 4
    for (int i = 0; i < 4; ++i) {
        local[i] = items[i].g;
    }

    for (int i = 0; i < 4; ++i) {
        const std::vector<double>& base_pos =
            global_nodes_coordinates.at(local[i]);
        double best = -1e308;
        int bestj = -1;
#pragma unroll 4
        for (size_t j = 0; j < others.size(); ++j) {
            if (used_opp[j] != 0) {
                continue;
            }
            std::vector<double> r = MathHelper::subtractVectors(
                global_nodes_coordinates.at(others[j]), base_pos);
            double val = MathHelper::calculateScalarProductOfVectors(r, n);
            if (val > best) {
                best = val;
                bestj = static_cast<int>(j);
            }
        }
        used_opp[bestj] = 1;
        local[4 + i] = others[bestj];
    }

    return local;
}

std::vector<std::vector<int>> Grid::applyTemplate(
    const std::vector<int>& local_mapping) {
    std::vector<std::vector<int>> elems;

    elems.reserve(5);
    elems.emplace_back();

#pragma unroll 4
    for (int k = 0; k < 5; ++k) {
        elems.back().push_back(local_mapping[PYRAMID_TEMPLATE.at(k)]);
    }

    for (int t = 0; t < 4; ++t) {
        elems.emplace_back();
#pragma unroll 4
        for (int k = 0; k < 4; ++k) {
            elems.back().push_back(
                local_mapping[TETRAHEDRON_TEMPLATE.at(t).at(k)]);
        }
    }

    return elems;
}

void Grid::ensurePositiveTets(
    std::vector<std::vector<int>>& elems,
    const std::map<int, std::vector<double>>& global_nodes_coordinates) {
#pragma unroll 4
    for (size_t e = 1; e < elems.size(); ++e) {
        if (elems[e].size() != 4) {
            continue;
        }
        const std::vector<double>& a = global_nodes_coordinates.at(elems[e][0]);
        const std::vector<double>& b = global_nodes_coordinates.at(elems[e][1]);
        const std::vector<double>& c = global_nodes_coordinates.at(elems[e][2]);
        const std::vector<double>& d = global_nodes_coordinates.at(elems[e][3]);
        double s6 = MathHelper::calculateOrientedTetrahedronVolume6(a, b, c, d);
        if (s6 < 0.0) {
            std::swap(elems[e][0], elems[e][1]);
        }
    }
}

void Grid::splitHexAnyFixedFace(
    const std::vector<int>& global_nodes_indicies,
    const std::map<int, std::vector<double>>& global_nodes_coordinates,
    const std::vector<int>& fixed_base_global_indicies) {
    auto local_mapping = buildLocalMappingForFixedFace(
        global_nodes_indicies, global_nodes_coordinates,
        fixed_base_global_indicies);

    auto elems = applyTemplate(local_mapping);

    ensurePositiveTets(elems, global_nodes_coordinates);

    pyramide_indicies.push_back(elems[0]);

#pragma unroll 4
    for (size_t i = 1; i < elems.size(); i++) {
        tetrahedron_indicies.push_back(elems[i]);
    }
}

void Grid::printElements() {
    std::cout << "Pyramids: \n";

    for (const auto& pyramide_index : pyramide_indicies) {
        std::cout << "\t[";
#pragma unroll 4
        for (size_t i = 0; i < pyramide_index.size(); ++i) {
            if (i != 0U) {
                std::cout << ", ";
            }
            std::cout << pyramide_index[i];
        }
        std::cout << "]\n";
    }

    std::cout << "Tetrahedrons: \n";

    for (const auto& tetrahedron_index : tetrahedron_indicies) {
        std::cout << "\t[";
#pragma unroll 4
        for (size_t i = 0; i < tetrahedron_index.size(); ++i) {
            if (i != 0U) {
                std::cout << ", ";
            }
            std::cout << tetrahedron_index[i];
        }
        std::cout << "]\n";
    }
}

void Grid::visualizeSplittedGrid(
    const std::map<int, std::vector<double>>& global_nodes_coordinates) {
    auto ugrid = VTKVisualizer::buildVTKGrid(
        global_nodes_coordinates, pyramide_indicies, tetrahedron_indicies);

    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0.9, 0.9, 0.9);

    auto pyr_filter = vtkSmartPointer<vtkThreshold>::New();
    pyr_filter->SetInputData(ugrid);
    pyr_filter->SetInputArrayToProcess(
        0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "elemType");
    pyr_filter->SetLowerThreshold(1);
    pyr_filter->SetUpperThreshold(1);
    pyr_filter->Update();

    auto pyr_surface = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    pyr_surface->SetInputConnection(pyr_filter->GetOutputPort());
    pyr_surface->Update();

    auto pyr_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    pyr_mapper->SetInputConnection(pyr_surface->GetOutputPort());

    auto pyr_actor = vtkSmartPointer<vtkActor>::New();
    pyr_actor->SetMapper(pyr_mapper);

    pyr_actor->GetProperty()->SetColor(1.0, 0.6, 0.2);
    pyr_actor->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
    pyr_actor->GetProperty()->SetOpacity(0.6);
    pyr_actor->GetProperty()->EdgeVisibilityOn();
    pyr_actor->GetProperty()->SetLineWidth(2.0);

    renderer->AddActor(pyr_actor);

    auto tet_filter = vtkSmartPointer<vtkThreshold>::New();
    tet_filter->SetInputData(ugrid);
    tet_filter->SetInputArrayToProcess(
        0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "elemType");
    tet_filter->SetLowerThreshold(2);
    tet_filter->SetUpperThreshold(2);
    tet_filter->Update();

    auto tet_surface = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    tet_surface->SetInputConnection(tet_filter->GetOutputPort());
    tet_surface->Update();

    auto tet_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    tet_mapper->SetInputConnection(tet_surface->GetOutputPort());

    auto tet_actor = vtkSmartPointer<vtkActor>::New();
    tet_actor->SetMapper(tet_mapper);
    tet_actor->GetProperty()->SetColor(0.2, 1.0, 0.2);
    tet_actor->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0);
    tet_actor->GetProperty()->SetOpacity(0.6);
    tet_actor->GetProperty()->EdgeVisibilityOn();
    tet_actor->GetProperty()->SetLineWidth(2.0);

    renderer->AddActor(tet_actor);

    auto ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    ren_win->AddRenderer(renderer);
    ren_win->SetSize(900, 700);
    ren_win->SetWindowName("Hex split visualization");

    auto iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);

    ren_win->Render();
    iren->Start();
}
