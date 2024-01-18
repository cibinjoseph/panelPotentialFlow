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

    return vindFil(p1, p2, p, self.rc)
end

function vindFil(p1::Vector{TF}, p2::Vector{TF}, p::Vector{TF}, rc::TF) where (TF<: Float64)
    r1 = p-p1
    r2 = p-p2
    r0 = r1-r2

    r1xr2 = cross(r1, r2)
    r1xr2abs2 = dot(r1xr2, r1xr2)

    r1Unit = r1/norm(r1)
    r2Unit = r2/norm(r2)

    vel = 0.0
    if r1xr2abs2 > eps2
        vel = (r1xr2*inv4pi*dot(r0, r1Unit-r2Unit)) / sqrt((rc*norm(r0))^4.0+r1xr2abs2^2.0)
    end
    return vel
end

struct Hshoe{TF}
    p::Matrix{TF}
    rc::TF
end

""" Constructor """
function Hshoe(p1, p2, p3, p4, rc)
    p = zeros(3, 4)
    p[:, 1] = p1
    p[:, 2] = p2
    p[:, 3] = p3
    p[:, 4] = p4
    return Hshoe(p, rc)
end

function vind(self::Hshoe, p)
    return (vindFil(self.p[:, 1], self.p[:, 2], p, self.rc) +
            vindFil(self.p[:, 2], self.p[:, 3], p, self.rc) +
            vindFil(self.p[:, 3], self.p[:, 4], p, self.rc))
end

struct Surface{TI, TF}
    mesh::TriMesh
    ele::Vector{Tri}
    nElem::TI
    gamma::Vector{TF}
    coeffP::Vector{TF}
    vel::Matrix{TF}
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
    for i = 1:nElem
        ele[i] = Tri(mesh.points[mesh.cells[i]][1],
                     mesh.points[mesh.cells[i]][2],
                     mesh.points[mesh.cells[i]][3], rc)
    end

    # Compute aic matrix
    aic = zeros(nElem, nElem)

    for j = 1:nElem
        for i = 1:nElem
            aic[i, j] = dot(vind(ele[j], ele[i].cp), ele[i].ncap)
        end
    end

    coeffP = zeros(nElem)
    vel = zeros(3, nElem)
    return Surface(mesh, ele, nElem, gamma, coeffP, vel, isClosed, aic)
end

function vind(self::Surface, p)
    vel = zeros(3)
    for i = 1:self.nElem
        vel += vind(self.ele[i], p) .* self.gamma[i]
    end
    return vel
end

function set!(RHS::Vector{Float64}, self::Surface, vinf::Vector{Float64}, hshoes::Vector{Hshoe}, gamma_hshoes::Vector{Float64})
    for i = 1:self.nElem
        vel_wake = vind(hshoes, gamma_hshoes, self.ele[i].cp)
        RHS[i] = -1 .* dot(vinf + vel_wake, self.ele[i].ncap)
    end
end

function vind(hshoes::Vector{Hshoe}, gamma_hshoes::Vector{Float64}, p::Vector{Float64})
    vel = zeros(3)
    for ihs = 1:length(hshoes)
        vel += vind(hshoes[ihs], p) * gamma_hshoes[ihs]
    end
    return vel
end

function writeMesh(s::Surface, filename::String)
    points, cells = getVTKElements(s.mesh)
    ncap = zeros(3, s.nElem)
    cp = zeros(3, s.nElem)
    id = 1:s.nElem
    vel = zeros(3, s.nElem)
    for i = 1:s.nElem
        ncap[:, i] = s.ele[i].ncap
        cp[:, i] = s.ele[i].cp
    end
    vtk_grid(filename, points, cells) do vtk
        vtk["gamma"] = s.gamma
        vtk["cell_id"] = id
        vtk["cp", VTKCellData()] = cp
        vtk["ncap", VTKCellData()] = ncap
        vtk["vel", VTKCellData()] = s.vel
        vtk["coeffP", VTKCellData()] = s.coeffP
    end
end

function writeWake(hshoes, gamma_hshoes, filename)
    n = length(hshoes)
    xyz = zeros(3, 2, n+1, 1)

    for i = 1:n
        xyz[:, 1, i, 1] .= hshoes[i].p[:, 2]
        xyz[:, 2, i, 1] .= hshoes[i].p[:, 1]
    end
    xyz[:, 1, n+1, 1] .= hshoes[n].p[:, 3]
    xyz[:, 2, n+1, 1] .= hshoes[n].p[:, 4]

    vtk_grid(filename, xyz) do vtk
        vtk["gamma"] = gamma_hshoes
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
    for i = 1:np
        if s.mesh.points[i][1] > max_x
            max_x = s.mesh.points[i][1]
            nodeList = [i]
        elseif abs(s.mesh.points[i][1] - max_x) < 1.0e-6
            push!(nodeList, i)
        end
    end

    # Identify cells that share these nodes
    cellsTE = zeros(Int, 2*(length(nodeList)-1))
    cellsU::Vector{Int} = []
    cellsL::Vector{Int} = []
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

    # Identify upper and lower cells based on ncap value
    zcap = [0, 0, 1.0]
    for i in cellsTE
        zcomp = dot(s.mesh.normals[i], zcap)
        if zcomp > 0.0
            push!(cellsU, i)
        else
            push!(cellsL, i)
        end
    end
    return nodeList, cellsTE, cellsU, cellsL
end

function set_hshoe_tail!(p1, p2, p3, p4; extend=3.0)
    p1 .= p2 .+ [extend*p2[1], 0.0, 0.0]
    p4 .= p3 .+ [extend*p3[1], 0.0, 0.0]
end

function get_hshoes(s::Surface; rc=1.0e-6)
    nodeList, cellsTE, cellsU, cellsL = getTE(s)
    hshoes = Vector{Hshoe}(undef, length(nodeList)-1)

    for i = 1:length(nodeList)-1
        p1 = zeros(3)
        p2 = s.mesh.points[nodeList[i]]
        p3 = s.mesh.points[nodeList[i+1]]
        p4 = zeros(3)
        set_hshoe_tail!(p1, p2, p3, p4)
        hshoes[i] = Hshoe(p1, p2, p3, p4, rc)
    end
    return hshoes, cellsU, cellsL
end

function set!(gamma_hshoes::Vector{Float64}, s::Surface, cellsU, cellsL)
    for i = 1:length(gamma_hshoes)
        gamma_hshoes[i] = s.gamma[cellsU[i]] - s.gamma[cellsL[i]]
    end
end

function setLocalVel!(s::Surface, hshoes::Vector{Hshoe}, gamma_hshoes::Vector{Float64}, vinf=Vector{Float64})
    for i = 1:s.nElem
        s.vel[:, i] = vinf + vind(s, s.ele[i].cp) + vind(hshoes, gamma_hshoes, s.ele[i].cp)
        s.coeffP[i] = norm(s.vel[:, i])^2
    end
    s.coeffP .= 1.0 .- s.coeffP ./ norm(vinf)^2
end

# Main program
body = Surface("geometry/airfoil.stl"; isClosed=false)
RHS = zeros(body.nElem)
vinf = [1.0, 0, 0]
hshoes, cellsU, cellsL = get_hshoes(body)
gamma_hshoes = zeros(length(hshoes))

# Initial iteration without wake
set!(RHS, body, vinf, hshoes, gamma_hshoes)
gamma = solve(body.aic, RHS; isClosed=body.isClosed)
assignGamma!(body, gamma)
set!(gamma_hshoes, body, cellsU, cellsL)
gamma_prev = gamma

# Iterate till wake and body strengths are converged
for iter = 1:50
    set!(RHS, body, vinf, hshoes, gamma_hshoes)
    gamma = solve(body.aic, RHS; isClosed=body.isClosed)
    assignGamma!(body, gamma)
    set!(gamma_hshoes, body, cellsU, cellsL)
    residual = norm(gamma - gamma_prev)
    println("Iter: $iter  Res: $residual")
    gamma_prev .= gamma

    if residual < 1e-6; break; end
end

setLocalVel!(body, hshoes, gamma_hshoes, vinf)
writeMesh(body, "out")
writeWake(hshoes, gamma_hshoes, "out_wake")
