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


pgfplotsx()
t, y = generate_signal()
plot(t, y);
for k = 1:3
  t, y = generate_signal()
  plot!(t, y);
end
plot!(show = true, legend=nothing, xlabel=latexify("ω"), ylabel=latexify("y(ω)"))
savefig("../TeX/PlotSources/less_terrible.pdf")
