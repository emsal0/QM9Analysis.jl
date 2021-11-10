using Flux
using Flux: @epochs
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
            CGConv((4, 1)),
            CGConv((4, 1)),
            CGConv((4, 1)),
            CGConv((4, 1)),
            node_feature,
            X -> sum(X, dims=2),
            Dense(4,10), Dense(10,10),
            Dense(10,4), Dense(4,1)
           )

loss(gf, y) = Flux.Losses.mse(model(gf), y)

NUM_TRAIN = convert(Int, floor(length(dat) * 4/5))

dat_train = [(row["molecule"], parse(Float64, row["R2"])) for row in dat[1:NUM_TRAIN]]
dat_val  = [(row["molecule"], parse(Float64, row["R2"])) for row in dat[NUM_TRAIN+1:end]]

@epochs 10 Flux.train!(loss, Flux.params(model), dat_train, ADAM(0.05))

losses_val = [loss(x,y) for (x,y) in dat_val]
avg_loss_val = 1/length(losses_val) * sum(losses_val)
@show max(losses_val)
@show avg_loss_val

@save "model1.bson" model
