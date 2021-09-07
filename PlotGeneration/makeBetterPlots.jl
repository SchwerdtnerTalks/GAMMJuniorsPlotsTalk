using LinearAlgebra, Plots, Statistics, SignalAnalysis, Latexify
plotly()
function generate_signal(;
    n=1000,
    t = range(0, 10, length=n)
  )
  f1(x) = cos(x)
  signal = f1.(t) + rand(PinkGaussian(length(t))) .* sin.(2t) .* t
  return t, signal
end

function generate_better_signal(;
    n=1000,
    t = range(0, 10, length=n),
    α=1.5
  )
  f1(x) = 2*sqrt(x)+cos(log(x+0.1))
  signal = f1.(t) + rand(PinkGaussian(length(t)))*α + rand(PinkGaussian(length(t))) .* α .* sin.(2t)
  return t, signal
end

struct SignalCollection{TSV}
  SV::TSV
  t
end

SV = Vector{Vector{Float64}}(undef, 0)
t, y = generate_signal()
for k = 1:1000
  t, y = generate_signal()
  push!(SV, y)
end
SC = SignalCollection(SV, Vector(t))

@recipe function f(r::SignalCollection; ε_max = 0.5)
  # set a default value for an attribute with `-->`
  xguide --> latexify("ω")
  yguide --> latexify("y(ω)")
  # add a series for an error band
  m, ep = mean_std(r)
  @series begin
    # force an argument with `:=`
    seriestype := :path
    # ignore series in legend and color cycling
    primary := false
    linecolor := nothing
    fillcolor := :lightgray
    fillalpha := 0.5
    fillrange := m .- ep
    # ensure no markers are shown for the error band
    markershape := :none
    # return series data
    r.t, m .+ ep
  end
  # get the seriescolor passed by the user
  c = get(plotattributes, :seriescolor, :auto)
  # highlight big errors, otherwise use the user-defined color
  seriescolor := ifelse.(ep .> ε_max, :red, c)
  # return data
  r.t, m
end

function mean_std(sc::SignalCollection)
  ep = Vector{Float64}(undef, length(sc.SV[1]))
  m = Vector{Float64}(undef, length(sc.SV[1]))
  v = Vector{Float64}(undef, length(sc.SV))
  for i = eachindex(ep)
    for j in eachindex(v)
      v[j] = sc.SV[j][i]
    end
    ep[i] = std(v)
    m[i] = mean(v)
  end
  return m, ep
end


SV = Vector{Vector{Float64}}(undef, 0)
t, y = generate_better_signal(α=0.9)
for k = 1:1000
  t, y = generate_better_signal(α=0.9)
  push!(SV, y)
end
SCnew = SignalCollection(SV, Vector(t))

pgfplotsx()
plot(SC, ε_max=50, legend=nothing)
savefig("../TeX/PlotSources/std1.pdf")

plot(SC, ε_max=5, legend=nothing)
savefig("../TeX/PlotSources/std2.pdf")

plot(SC, ε_max=5, label="Miller et al.", show=false);
plot!(SCnew, ε_max=5, seriescolor=:green, label = "our method", legend = :topright)
savefig("../TeX/PlotSources/std3.pdf")

plot(SC, ε_max=5, label="Miller et al.")
plot!(SCnew, ε_max=5, seriescolor=:green, label = "our method", legend = :bottomleft)
savefig("../TeX/PlotSources/std4.pdf")
