using Flux
using GeometricFlux
using Serialization
using BSON: @save

DATAFILE="data.blob"
dat = open(DATAFILE) do io
    deserialize(io)
end

# @show dat
#
model = Flux.Chain(
            CGConv((4, 1)),
            GINConv(Flux.Chain(Dense(4,5), Dense(5,5), Dense(5, 3))),
            node_feature,
            X -> sum(X, dims=2),
            Dense(3,4), Dense(4,4), Dense(4,1)
           )

loss(gf, y) = Flux.Losses.mse(model(gf), y)

dat_train = [(row["molecule"], parse(Float64, row["R2"])) for row in dat]

Flux.train!(loss, Flux.params(model), dat_train, ADAM(0.01))

@save "model1.bson" model
