mutable struct VertexList{Model}
    vertices::Array{Int}
    v_first::Matrix{Int}
    v_last::Matrix{Int}
end

function VertexList{Model}(site_count::Integer) where {Model}
    return VertexList(
        zeros(Int, 0, 0, 0),
        zeros(Int, 2, site_count),
        zeros(Int, 2, site_count),
    )
end

function make_vertex_list!(
    vl::VertexList{Model},
    operators::AbstractVector{<:OperCode},
    data::SSEData,
) where {Model}
    if size(vl.vertices, 2) != length(operators)
        vl.vertices = Array{Int,3}(undef, 2, leg_count(Model), length(operators))
    end

    vl.vertices .= -1
    vl.v_first .= -1
    vl.v_last .= -1

    for (p, op) in enumerate(operators)
        if isidentity(op)
            continue
        end

        b = data.bonds(get_bond(op))
        for s = 1:leg_count(Model)÷2
            (s1, p1) = vl.v_last[:, b.sites[s]]
            if p1 != -1
                vl.vertices[:, s1, p1] .= (s, p)
                vl.vertices[:, s, p] .= (s1, p1)
            else
                vl.v_first[:, b.sites[s]] .= (s, p)
            end
            vl.v_last[:, b.sites[s]] .= (leg_count(Model) ÷ 2 + s, p)
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
