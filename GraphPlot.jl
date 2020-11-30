using Plots

using DelimitedFiles

function f(A,x)
	sum=0
	for i=1:size(A,1)
	sum+=A[i,1]*(x-x_a)^(i-1)
	end
	return sum
end

for i=1:2796
	A = readdlm("SL/"*list[i,1])
	B = readdlm("Poly/"*list[i,1])
	C = readdlm("res/"*list[i,1])
	p = scatter(C[:,1],C[:,2],legend=false)
	g(x) = f(B,x)
	plot!(g,[0.05:0.01:1])
	lower = upper = 0
	for j=1:size(A,1)
		lower = upper
		upper = A[j,1]
		h(x) = A[j,2]*x+A[j,3]
		interval = (upper-lower)/10
		plot!(h,[lower,interval,upper])
	end
	savefig("png/"*list[i,1])
end



for j=5:9
	zeros000="0."
	for k=1:j-1
		zeros000 = zeros000*"0"
	end
	for i=1:2796
		A=list[i,1]
		if A[length(A)-10:length(A)-8]=="e-"*string(j)
			mv(A,zeros000*A[1]*A[3:length(A)-11]*".png")
		end
	end
end