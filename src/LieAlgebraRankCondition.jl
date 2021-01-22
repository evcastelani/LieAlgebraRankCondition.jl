
module LieAlgebraRankCondition

	using LinearAlgebra

	export lieprod,isli,rankcondition,ad, OutputSet

	import Base.show

	"""
		OutputSet
	
	It is an internal type used to print the output of `rankcondition` function. Once the `rankcondition` function has been used, just access the fields according to the example below.
	
	# Examples
	```
	julia-repl
	julia> s = rankcondition(A,B)
	julia> s.LI # returns true or false
	julia> s.elements # returns a sequence of arrays (L.I arrays)
	julia> s.path # show the sequence of Lie products used to get the solution `s`.
	```
	"""	
	struct OutputSet
		path :: Vector{String}
		LI :: Bool
		elements :: Vector
	end

	"""
		lieprod
	
	An auxiliary function to compute the lie product of two matrix: AB-BA
	
	# Examples
	```
	julia-repl
	julia> A = rand(4,4)
	julia> B = rand(4,4)
	julia> C = lieprod(A,b)
	```
	returns the lie product of A and B and storage in C.
	"""
	function lieprod(A::Array{Float64},B::Array{Float64})
		return A*B-B*A
	end

	"""
		isli
	
	An auxiliary function to test if a set of matrix is an L. I set.
	
	# Examples
	```
	julia-repl
	julia> C=[A,B]
	julia> t = isli(C)
	```
	returns true or false and storage the answer in t.
	"""
	function isli(C::Vector)
		D = zeros(length(C),16)
		for k=1:length(C)
			D[k,:]=C[k]
		end
		if rank(D)==length(C)
			return true
		else
			return false
		end
	end
	"""
		rankcondition
	
	It is the main function. See the example to run. 
	
	# Example
	```
	julia-repl
	julia> A = [1 1 0 0; -1.0 1 0 0;0 0 -1 0.5;0 0 -0.5 -1]
	julia> B = [2.0 0.0 0.0 0.0; 0.0 -1.5 -0.1 0.0; 0.0 0.1 -1.5 0.0; 0.0 0.0 0.0 1.0]
	julia> s = rankcondition(A,B)
	```
	returns s like an OutputSet with three fields. This OutputSet type is used to store L.I. matrices, bool values and a convenient path to reproduce the same solution. There is an optional argument used to define the number of L.I. elements. 
	# Example
	```
	julia-repl
	julia> A = [1 1 0 0; -1.0 1 0 0;0 0 -1 0.5;0 0 -0.5 -1]
	julia> B = [ 2.0   0.0   0.0  0.0; 0.0  -1.5  -0.1  0.0; 0.0   0.1  -1.5  0.0; 0.0   0.0   0.0  1.0]
	julia> s = rankcondition(A,B,dim=15)
	```
	By default `dim` parameter can be omitted. You can change this parameter if it is convenient for you. 
	"""
	function rankcondition(A,B,dim=15)
		partial = [A,B]
		path = ["A","B"]
		search_index = ones(Int,dim)
		if isli(partial)==true
			push!(partial,lieprod(partial[1],partial[2]))
			push!(path,"[$(path[1]),$(path[2])]")
			k = 3 
			while 3<=k<=dim
				next_step = false
				search_partial_solution = true
				while search_partial_solution == true
					if isli(partial)==true
						push!(partial,lieprod(partial[search_index[k]],partial[k]))
						push!(path,"[$(path[search_index[k]]),$(path[k])]")
						next_step = true
						search_partial_solution = false
					elseif search_index[k]<=k-1
						search_index[k] = search_index[k]+1
						partial[k] = lieprod(partial[search_index[k]],partial[k-1])
						path[k] = "[$(path[search_index[k]]),$(path[k-1])]"
						search_partial_solution = false
					end
				end
				if next_step == false
					deleteat!(partial,k)
					deleteat!(path,k)
					search_index[k] = 1
					k = k-1
					search_index[k]=search_index[k]+1
					partial[k] = lieprod(partial[search_index[k]],partial[k-1])
					path[k] = "[$(path[search_index[k]]),$(path[k-1])]"
				else
					k = k+1
				end
			end
			if k<=3
				println(" 🚑 No L.I sets were found!")
				return OutputSet([""],false,[])
			end
			if isli(partial[1:dim]) == true
				println(" 🎈 A set of L.I arrays with $(dim) elements was found.")
				return OutputSet(path[1:dim],true,partial[1:dim])
			end
		else 
			println(" 🚑 No L.I sets with $(dim) elements were found!")
			return OutputSet([""],false,[])

		end
	end

	"""
		ad(n,A,B)
	
		This functins is used to compute [A[A,...,[A,B]]
	
	# Examples
	```
	julia-repl
	julia> ad(2,A,B)
	```
	returns lieprod(A,lieprod(A,B))
	"""
	function ad(n,M1,M2)
		M = lieprod(M1,M2)
		for k=2:n
			M = lieprod(M1,M)
		end
		return M
	end

	#####################
	# print output	    #			                                                                    #
	#####################

	function Base.show(io::IO,ots::OutputSet)
		if ots.LI == true
			print(io,"\n 👉 Here are the set and the path:  \n")
			for k = 1:length(ots.elements)
				print(io, "\n ↪  Element $(k):")
				print(io, " $(ots.path[k]) = \n")
				display(ots.elements[k])
			end
		else
			print(io, "\n ")
		end
	end
end # module
