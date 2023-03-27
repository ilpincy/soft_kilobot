using CSV, DataFrames, Plots


rodi = CSV.read("rodi.txt", DataFrame; delim=" ");
polar = CSV.read("polar.txt", DataFrame; delim=" ");
diff = DataFrame(
  p_x = rodi[!,:p_x],
  p_y = rodi[!,:p_y],
  f = rodi[!,:f] .- polar[!,:f]
);
xs = unique(diff[!,:p_x])
ys = unique(diff[!,:p_y])
data = zeros(length(xs), length(ys))
xsi = Dict(x => i for (i,x) in enumerate(xs))
ysi = Dict(y => i for (i,y) in enumerate(ys))
for row in 1:size(diff,1)
  p_x = diff[row, :p_x];
  p_y = diff[row, :p_y];
  f = diff[row, :f];
  data[xsi[p_x], ysi[p_y]] = f;
end

maxdiff = maximum(abs.(diff[!,:f]));
println("The maximum difference is ", maxdiff);

gr();
plt = heatmap(xs,ys,data,size=(600,600));
savefig(plt, "comparison.pdf");
display(plt);
