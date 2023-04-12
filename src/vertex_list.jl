mutable struct VertexList{NSites}
    vertices::Array{Int}
    v_first::Matrix{Int}
    v_last::Matrix{Int}
end

function VertexList{NSites}(site_count::Integer) where {NSites}
    return VertexList{NSites}(
        zeros(Int, 0, 0, 0),
        zeros(Int, 2, site_count),
        zeros(Int, 2, site_count),
    )
end

function make_vertex_list!(
    vl::VertexList{NSites},
    operators::AbstractVector{<:OperCode},
    bonds::AbstractVector{<:Bond{NSites}},
) where {NSites}
    if size(vl.vertices, 2) != length(operators)
        vl.vertices = fill(-1, 2, 2 * NSites, length(operators))
    end

    vl.v_first .= -1
    vl.v_last .= -1

    for (p, op) in enumerate(operators)
        if isidentity(op)
            continue
        end

        b = bonds[get_bond(op)]
        for s = 1:NSites
            (s1, p1) = vl.v_last[:, b.sites[s]]
            if p1 != -1
                vl.vertices[:, s1, p1] .= (s, p)
                vl.vertices[:, s, p] .= (s1, p1)
            else
                vl.v_first[:, b.sites[s]] .= (s, p)
            end
            vl.v_last[:, b.sites[s]] .= (NSites + s, p)
        end
    end

    for i = 1:size(vl.v_first, 2)
        if vl.v_first[1, i] != -1
            vl.vertices[:, vl.v_first[:, i]...] .= vl.v_last[:, i]
            vl.vertices[:, vl.v_last[:, i]...] .= vl.v_first[:, i]
        end
    end

    return nothing
end
