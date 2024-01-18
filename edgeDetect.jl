using VSPGeom
using VSPGeom.WriteVTK
using LinearAlgebra

inv4pi = 0.07957747155
eps2 = eps()^2

""" Class for handling triangular vortex elements """
struct Tri{TF}
    p::Matrix{TF}
    ncap::Vector{TF}
    rc::TF
    cp::Vector{TF}
end

""" Constructor """
function Tri(p, rc)
    ncap = cross(p[:, 2]-p[:, 1], p[:, 3]-p[:, 2])
    ncap .= ncap/norm(ncap)

    cp = zeros(3)
    cp .= sum(p; dims=2) ./ 3
    return Tri(p, ncap, rc, cp)
end

function Tri(p1::Vector{TF}, p2::Vector{TF}, p3::Vector{TF}, rc::TF) where (TF<:Float64)
    p = zeros(3, 3)
    p[:, 1] = p1
    p[:, 2] = p2
    p[:, 3] = p3

    return Tri(p, rc)
end

struct Surface
    mesh::TriMesh
    ele::Vector{Tri}
    nElem::Int64
    isClosed::Bool
    isTrailing::Vector{Int64}
end

""" Constructor """
function Surface(filename::String; rc=1e-6, isClosed=true, tol=1e-6)
    geom = readSTL(filename)
    mesh = geom[1]
    nElem = length(mesh.cells)
    gamma = zeros(nElem)

    ele = Vector{Tri}(undef, nElem)
    for i = 1:nElem
        ele[i] = Tri(mesh.points[mesh.cells[i]][1],
                     mesh.points[mesh.cells[i]][2],
                     mesh.points[mesh.cells[i]][3], rc)
    end
    isTrailing = zeros(nElem)
    return Surface(mesh, ele, nElem, isClosed, isTrailing)
end

function writeMesh(s::Surface, filename::String)
    points, cells = getVTKElements(s.mesh)
    ncap = zeros(3, s.nElem)
    cp = zeros(3, s.nElem)
    id = 1:s.nElem
    vel = zeros(3, s.nElem)
    isTrailing = zeros(s.nElem)

    for i = 1:s.nElem
        ncap[:, i] = s.ele[i].ncap
        cp[:, i] = s.ele[i].cp
        isTrailing[i] = s.isTrailing[i]
    end

    vtk_grid(filename, points, cells) do vtk
        vtk["cell_id"] = id
        vtk["cp", VTKCellData()] = cp
        vtk["ncap", VTKCellData()] = ncap
        vtk["trailing"] = isTrailing
    end
end

function getNeighborCells(s::Surface, icell)
    neighborCells = zeros(Int, 3)
    vtxs = s.mesh.cells[icell]
    for i in 1:s.nElem
        if vtxs[1] in s.mesh.cells[i] && vtxs[2] in s.mesh.cells[i]
            neighborCells[1] = i
        elseif vtxs[2] in s.mesh.cells[i] && vtxs[3] in s.mesh.cells[i]
            neighborCells[2] = i
        elseif vtxs[3] in s.mesh.cells[i] && vtxs[1] in s.mesh.cells[i]
            neighborCells[3] = i
        end
    end
    return neighborCells
end

function angleBetween_deg(a, b)
    return atand(norm(cross(a, b)), dot(a, b))
end

function getTrailingNodes(s::Surface; Vinf_dir=[1.0, 0.0, 0.0], tol_deg=20)
    cellidx::Vector{Int64} = []
    edgeNodes::Vector{Int64} = []
    edgeStart = [1, 2, 3]
    edgeEnd = [2, 3, 1]

    for i = 1:s.nElem
        # if i == 288
        # For each cell, if the normal does not
        # have a component along Vinf, ignore it.
        ncap = s.ele[i].ncap
        if dot(ncap, Vinf_dir) >= 0.0
            neighborCells = getNeighborCells(s, i)
            # Find the cell that shares an edge
            for iedge = 1:3
                nearCell = neighborCells[iedge]
                # Check the angle between normal vectors
                # Only check cells with index larger than the current one
                # as the others would already have been checked
                if nearCell > i && angleBetween_deg(s.ele[nearCell].ncap, ncap) > tol_deg
                    push!(edgeNodes, edgeStart[iedge])
                    push!(edgeNodes, edgeEnd[iedge])
                    push!(cellidx, i)
                end
            end
        end
        # end
    end
    return unique(cellidx), unique(edgeNodes)
end

# Main program
body = Surface("airfoil.stl"; isClosed=false)
cellidx, edgeNodes = getTrailingNodes(body)
body.isTrailing[cellidx] .= 1
writeMesh(body, "edge")
