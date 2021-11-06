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

for i in 1:length(dat)
    println(model(dat[i]["molecule"]), " ", dat[i]["R2"])
end
