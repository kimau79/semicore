Date last updated: 7 Jan 2020

Command to run the code from Unix:
"julia main.jl (input1) (input2) (input3)"
(input1) = Thickness of first shell (see description of "firstShellThickness" from main.jl)
(input2) = shell thickness factor (see description of "shellThicknessFactor" from main.jl)
(input3) = order of polynomial (see "orderOfpolynomial" from main.jl)

Example:
julia main.jl 1e-7 1.003 14

NOTE:
In the programming of main.jl,
(input1) corresponds to "parse(Float64,ARGS[1])"
(input2) corresponds to "parse(Float64,ARGS[2])"
(input3) corresponds to "parse(Int64,ARGS[3])"
You might change the above parts for easier command input at you wish

All intermediary outputs for my own checking (e.g. outputs of WeightFactor) has been REMOVED.
