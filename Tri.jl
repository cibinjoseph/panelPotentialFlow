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

""" Induced velocity """
function vind(self::Tri, p)
    return vindFil(self, 1, p) + vindFil(self, 2, p) + vindFil(self, 3, p)
end


""" Induced velocity by filament """
function vindFil(self::Tri, ifil::Int, p::Vector{Float64})
    if ifil < 3
        p1 = self.p[:, ifil]
        p2 = self.p[:, ifil+1]
    else
        p1 = self.p[:, 3]
        p2 = self.p[:, 1]
    end

    r1 = p-p1
    r2 = p-p2
    r0 = r1-r2

    r1xr2 = cross(r1, r2)
    r1xr2abs2 = dot(r1xr2, r1xr2)

    r1Unit = r1/norm(r1)
    r2Unit = r2/norm(r2)

    vel = 0.0
    if r1xr2abs2 > eps2
        vel = (r1xr2*inv4pi*dot(r0, r1Unit-r2Unit)) / sqrt((self.rc*norm(r0))^4.0+r1xr2abs2^2.0)
    end
    return vel
end

struct Surface{TI, TF}
    mesh::TriMesh
    ele::Vector{Tri}
    nElem::TI
    gamma::Vector{TF}
    isClosed::Bool
    aic::Matrix{TF}
end

""" Constructor """
function Surface(filename::String; rc=1e-6, isClosed=true, tol=1e-6)
    geom = readSTL(filename)
    mesh = geom[1]
    nElem = length(mesh.cells)
    gamma = zeros(nElem)

    ele = Vector{Tri}(undef, nElem)
    for i in 1:nElem
        ele[i] = Tri(mesh.points[mesh.cells[i]][1],
                     mesh.points[mesh.cells[i]][2],
                     mesh.points[mesh.cells[i]][3], rc)
    end

    # Compute aic matrix
    aic = zeros(nElem, nElem)

    for j in 1:nElem
        for i in 1:nElem
            aic[i, j] = dot(vind(ele[j], ele[i].cp), ele[i].ncap)
        end
    end

    return Surface(mesh, ele, nElem, gamma, isClosed, aic)
end

function vind(self::Surface, p)
    vel = 0.0
    for i in range(1:self.nElem)
        vel += vind(self.ele[i], p) .* self.gamma[i]
    end
    return vel
end

function getRHS(self::Surface, vinf)
    RHS = zeros(self.nElem)

    for i in 1:self.nElem
        RHS[i] = -1 .* dot(vinf, self.ele[i].ncap)
    end

    return RHS
end

function writeMesh(s::Surface, filename::String; vinf=zeros(3))
    points, cells = getVTKElements(s.mesh)
    ncap = zeros(3, s.nElem)
    cp = zeros(3, s.nElem)
    id = 1:s.nElem
    for i = 1:s.nElem
        ncap[:, i] = s.ele[i].ncap
        cp[:, i] = s.ele[i].cp
    end
    vtk_grid(filename, points, cells) do vtk
        vtk["gamma"] = s.gamma
        vtk["cell_id"] = id
        vtk["cp", VTKCellData()] = cp
        vtk["ncap", VTKCellData()] = ncap
        vtk["vel", VTKCellData()] = repeat(vinf, outer=(1, s.nElem))
    end
end

function solve(aic, RHS; isClosed=true)
    if isClosed
        return aic[2:end, 2:end] \ RHS[2:end]
    else
        return aic \ RHS
    end
end

function assignGamma!(s::Surface, gamma)
    idxStart = 1
    if length(gamma) == s.nElem-1
        s.gamma[1] = 0.0
        idxStart = 2
    end
    s.gamma[idxStart:end] .= gamma
end

function getTE(s::Surface)
    # Use max_x based search - not valid for swept geometry
    np = length(s.mesh.points)
    # Find max x coordinate
    max_x = s.mesh.points[1][1]
    nodeList::Vector{Int} = []
    for i in 1:np
        if s.mesh.points[i][1] > max_x
            max_x = s.mesh.points[i][1]
            nodeList = [i]
        elseif abs(s.mesh.points[i][1] - max_x) < 1.0e-6
            push!(nodeList, i)
        end
    end

    # Identify cells that share these nodes
    cellsTE = zeros(Int, 2*(length(nodeList)-1))
    cellsU = zeros(Int, length(nodeList)-1)
    cellsL = zeros(Int, length(nodeList)-1)
    j = 1
    for inode=1:length(nodeList)-1
        for i=1:s.nElem
            if nodeList[inode] in s.mesh.cells[i]
                if nodeList[inode+1] in s.mesh.cells[i]
                    cellsTE[j] = i
                    j += 1
                end
            end
        end
    end

    # Identify upper and lower cells based on z value
    for i in 1:length(cellsTE)
        crds = s.mesh.points[body.mesh.cells[i]]
        zs = [crds[1][3], crds[2][3], crds[3][3]]
        @show zmax = max(zs...)
        @show zmin = min(zs...)
    end
    return nodeList, cellsTE, cellsU, cellsL
end

function vind_chordwisewake(s::Surface)
    nodeList = getTE(s)
    # return vel
end

# Main program
body = Surface("airfoil.stl")
vinf = [1,0,0]
RHS = getRHS(body, vinf)
println("Computing solution ...")
gamma = solve(body.aic, RHS; isClosed=body.isClosed)
assignGamma!(body, gamma)
nodeList, cellsTE, cellsU, cellsL = getTE(body)
writeMesh(body, "out"; vinf=vinf)
