using VSPGeom
using VSPGeom.WriteVTK
using LinearAlgebra

struct Surface
    mesh::TriMesh
    nElem::Int64
    isTrailing::Vector{Int64}
end

""" Constructor """
function Surface(filename::String; tol=1e-6)
    geom = readSTL(filename)
    mesh = geom[1]
    nElem = length(mesh.cells)
    gamma = zeros(nElem)
    isTrailing = zeros(nElem)
    return Surface(mesh, nElem, isTrailing)
end

""" Writes out mesh to Paraview """
function writeMesh(s::Surface, filename::String)
    points, cells = getVTKElements(s.mesh)
    ncap = zeros(3, s.nElem)
    cp = zeros(3, s.nElem)
    id = 1:s.nElem
    vel = zeros(3, s.nElem)
    isTrailing = zeros(s.nElem)

    for i = 1:s.nElem
        ncap[:, i] = s.mesh.normals[i]
    end

    vtk_grid(filename, points, cells) do vtk
        vtk["cell_id", VTKCellData()] = id
        vtk["ncap", VTKCellData()] = ncap
        vtk["trailing", VTKCellData()] = s.isTrailing
    end
end

""" Obtain the 3 neighboring cells. 0 if there's no cell on that edge. """
function getNeighborCells(s::Surface, icell)
    neighborCells = zeros(Int, 3)
    vtxs = s.mesh.cells[icell]

    # Check if vertex is common to other cells' vertices
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

""" Angle between two vectors in degrees """
angleBetween_deg(a, b) = atand(norm(cross(a, b)), dot(a, b))

function getTrailingNodes(s::Surface; Vinf_dir=[1.0, 0.0, 0.0], edgeAngleTol_deg=40)
    cellidx::Vector{Int64} = []
    edgeNodes::Vector{Int64} = []
    edgeStart = [1, 2, 3]
    edgeEnd = [2, 3, 1]
    normalAngleTol_deg = 180.0-edgeAngleTol_deg

    for i = 1:s.nElem
        # For each cell, if the normal does not
        # have a component along Vinf, ignore it.
        ncap = s.mesh.normals[i]
        if dot(ncap, Vinf_dir) >= 0.0
            neighborCells = getNeighborCells(s, i)
            # Find the cell that shares an edge
            for iedge = 1:3
                nearCell = neighborCells[iedge]
                # Check the angle between normal vectors
                # Only check cells with index larger than the current one
                # as the others would already have been checked
                if nearCell > i && angleBetween_deg(s.mesh.normals[nearCell], ncap) > normalAngleTol_deg
                    push!(edgeNodes, edgeStart[iedge])
                    push!(edgeNodes, edgeEnd[iedge])
                    push!(cellidx, i)
                end
            end
        end
    end
    return unique(cellidx), unique(edgeNodes)
end

# Main program
body = Surface("geometry/aircraft.stl")
cellidx, edgeNodes = getTrailingNodes(body)
body.isTrailing[cellidx] .= 1
writeMesh(body, "edge")
