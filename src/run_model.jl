using BSON: @load
using Flux
using GeometricFlux
using Serialization

MODELFILE="model1.bson"
DATAFILE="data.blob"

@load MODELFILE model
dat = open(DATAFILE) do io
    deserialize(io)
end

loss(gf, y) = Flux.Losses.mse(model(gf), y)
for i in 1:length(dat)
    @show loss(dat[i]["molecule"], parse(Float64, dat[i]["R2"]))
end
