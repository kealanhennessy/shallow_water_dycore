module operators

using ..fields
using ..grid

function ∂λ!(df::Field, f::Field; g::Grid)
    @inbounds for i ∈ g.nx, j ∈ g.ny
        df.data[i,j] = f.data[i,j+1] - f.data[i,j] / g.dλ
    end
end

function ∂φ!()

end

function div!()

end

function grad!()

end

end