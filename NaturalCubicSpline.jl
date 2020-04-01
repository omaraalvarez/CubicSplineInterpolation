using Plots; gr();
using Plots.PlotMeasures
using LinearAlgebra;
using CSV;
using LaTeXStrings;
using ArgParse;

function main(File::String, Delimeter::String = ",")
    Data = CSV.read("$(pwd())/$(File)", delim = "$Delimeter", header = false, copycols = true);
    Data = hcat(Data.Column1, Data.Column2);
    issorted(Data[:, 1]) == false ? Data = Data[sortperm(Data[:, 1]), :] : nothing
    x, y = Data[:, 1], Data[:, 2];
    p = scatter((x, y), xlabel = L"x", ylabel = L"f(x)", color = "red", markersize = 6, legend = false, left_margin = 5mm, bottom_margin = 5mm, guidefontsize = 18, tickfontsize = 14, size = [1280, 720], dpi = 300)
    Interpolation!(x, y, File, p)
end

function Parse_Commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "File"
            help = "Name (and subroute) of the x, y file.";
            arg_type = String;
            required = true;
        "--Delimeter", "-D"
            help = "Data Delimiter (',' '\\t' ' ' etc.)";
            arg_type = String;
            required = false;
            default = ","
    end
    return parse_args(s)
end

function Interpolation!(x::Array{Float64}, y::Array{Float64}, File::String, p::Plots.Plot{Plots.GRBackend}, n::Int64 = 100)
    N, n = length(x), n*length(x);
    k = zeros(Float64, N);
    A = zeros(Float64, N-2, N-2);
    B = zeros(Float64, N-2);
    for i = 2:N - 1
        B[i - 1] = 6( (y[i-1] - y[i])/(x[i-1] - x[i]) - (y[i] - y[i+1])/(x[i] - x[i+1]) );
        if i != N-1
            A[i-1, i-1] = 2 * (x[i-1] - x[i+1]);
            A[i-1, i] = A[i, i-1] = x[i] - x[i+1];
        else
            A[i-1, i-1] = 2 * (x[i-1] - x[i+1]);
        end
    end
    SymTridiagonal(A);
    k_aux = inv(A) * B;
    for i = 1:N-2
        k[i+1] = k_aux[i];
    end
    E, i, r = x[1] + abs( (x[end] - x[1]) / n ), 1, 1;
    x_Interpolation, y_Interpolation = zeros(Float64, n), zeros(Float64, n);
    while E <= x[end]
        E >= x[i+1] ? i += 1 : nothing
        x_Interpolation[r] = E;
        α = (k[i] / 6) * (((E - x[i+1])^3) / (x[i] - x[i+1]) - (E - x[i+1]) * (x[i] - x[i+1]));
        β = (k[i+1] / 6)*(((E - x[i])^3) / (x[i] - x[i+1]) - (E - x[i]) * (x[i] - x[i+1]))
        γ = (y[i] * (E - x[i+1]) - y[i+1] * (E - x[i])) / (x[i] - x[i+1]);
        y_Interpolation[r] = α - β + γ;
        E += abs( (x[end] - x[1]) / n );
        r += 1;
    end
    if x_Interpolation[end] == 0.0
        resize!(x_Interpolation, n - 1)
        resize!(y_Interpolation, n - 1)
    end
    Output = open("$(pwd())/Interpolation_$(File[1:end-4]).dat", "w+")
    for i = 1:length(x_Interpolation)
        println(Output, "$(x_Interpolation[i])\t$(y_Interpolation[i])")
    end
    close(Output)

    plot!(p, (x_Interpolation, y_Interpolation), color = :black, width = 2)
    savefig(p, "$(pwd())/Interpolation_$(File[1:end-4]).png")
end

Args = Parse_Commandline()
println("Parsed args:")
for (arg,val) in Args
    println("  $arg  =>  $val")
end
@time main(Args["File"], Args["Delimeter"])
