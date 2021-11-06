using Flux
using GeometricFlux
using Serialization

DATAFILE="data.blob"
dat = open(DATAFILE) do io
    deserialize(io)
end

# @show dat

model = Flux.Chain(
            CGConv((4, 1)),
            GINConv(Flux.Chain(Dense(4,5), Dense(5, 1))),
           )

loss(gf, y) = Flux.Losses.mse(model(gf), y)

for row in dat
    molgf = row["molecule"]
    y = parse(Float64, row["R2"])
    @show loss(molgf, y)
end

