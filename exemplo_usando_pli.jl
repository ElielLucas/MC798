using JuMP
using Gurobi

function Exemplo_PLI(T::TriggerArcTSP, MaxTime::IntType)
	 global verbose_mode = true
	 global maximum_number_of_nodes_to_print
	(InArcs,OutArcs) = Build_In_Out_Arcs(T)
	if (verbose_mode)
	   Print_Incident_Arcs(T,InArcs,"InArc")
	   Print_Incident_Arcs(T,OutArcs,"OutArc")
	end

	if (verbose_mode)
   	   PrintHeader("Informacoes da execucao do solver")
	end
	TSP_model = Model(Gurobi.Optimizer)
	if (!verbose_mode)
	   set_optimizer_attribute(TSP_model, "OutputFlag", 0)

	   set_silent(TSP_model)
	end
	
	set_time_limit_sec(TSP_model,MaxTime)
	
	StartingTime = time()

	# create integer variables to for node sequence order
   	@variable(TSP_model, 1 <= y[1:T.NNodes] <= T.NNodes, Int)
	for v in 1:T.NNodes
	    set_name(y[v], "y_"*string(v))
	end

	# Node
	set_upper_bound(y[1],1) # y[1] fica fixado em 1
	# create binary variables for the arcs
   	@variable(TSP_model, x[i=1:T.NArcs], Bin)
	for a in 1:T.NArcs
	    set_name(x[a], "x_"*string(T.Arc[a].u)*"_"*string(T.Arc[a].v) )
	end

   	# In degree constraint:
	# Total sum of x_a, with arcs entering a node, is equal to 1.
   	for v in 1:T.NNodes
	    expr=0 # Start an empty expression
      	    for a in InArcs[v]
	    	 expr += x[a]
	    end
      	    @constraint(TSP_model, expr==1)
      	end

   	# Out degree constraint:
	# Total sum of x_a, with out leaving a node, is equal to 1.
   	for u in 1:T.NNodes
	    expr=0 # Start an empty expression
      	    for a in OutArcs[u]
	    	 expr += x[a]
	    end
      	    @constraint(TSP_model, expr==1)
      	end
   
	# We do not consider arcs entering node 1.
	for a in 1:T.NArcs
	    u = T.Arc[a].u
	    v = T.Arc[a].v
	    set_name(x[a], "x_"*string(u)*"_"*string(v) )
	    if (v!=1) # is an arc entering node 1.
               @constraint(TSP_model,y[u] + x[a] <= y[v] + (T.NNodes-1)*(1-x[a]))
            end
      	end
   	@objective(TSP_model, Min, sum(T.Arc[a].cost*x[a] for a in 1:T.NArcs))

	# for a in 1:T.NArcs
	#     if ((4,1) == (T.Arc[a].u,T.Arc[a].v))
	#        set_lower_bound(x[a],1) # Diz que vamos usar o arco (4,1)
	#     end
	# end

   	optimize!(TSP_model)
        if (termination_status(TSP_model) != MOI.OPTIMAL)
	   if (termination_status(TSP_model) == MOI.TIME_LIMIT)
	      error("It was not possible to obtain optimal solution, due to time limit.")
	   else 
              println(termination_status(TSP_model))
              error("It was not possible to obtain optimal solution.")
	   end
        end
	
	if (verbose_mode)
	   println("Obtained Optimal Solution")
   	   primal_status(TSP_model)
	   
	   # --------------------------------------------------
	   if (T.NNodes > maximum_number_of_nodes_to_print)
	      println("Verbose mode: Formulation is printed only for graphs up to ",maximum_number_of_nodes_to_print," nodes.")
	   else
	      PrintHeader("Formulacao MTZ para o TSP")
   	      print(TSP_model)
   	      # println("\nO modelo foi gravado no arquivo TSP_model.lp")
   	      # write_to_file(TSP_model, "TSP_model.lp")
	   end

	   # --------------------------------------------------
	   PrintHeader("Statistics")
	   
	   ExecutionTime = time() - StartingTime
	   println("Tempo maximo de execucao: ",MaxTime," seconds.")
	   println("Tempo de execucao usado:  ",ExecutionTime," seconds.")
	   
	   NumberOfNodes = MOI.get(TSP_model, MOI.NodeCount())
	   println("Number of nodes in the Branch and Bound tree: ",NumberOfNodes)
	   
	   obj_val = objective_value(TSP_model)
	   println("Value of the solution: ",obj_val)

	   ub_mtz::Vector{IntType} = Vector{IntType}(undef,T.NNodes)
	   
	   for u in 1:T.NNodes
      	       ub_mtz[u] = -1  # indice do arco da solucao que sai do vertice u
   	   end
   	   for a in 1:T.NArcs
               x_a = value(x[a])
	       if (x_a>0.0001)
	       	  if ub_mtz[T.Arc[a].u] != -1
	       	     error("Solucao define pelo menos dois arcos saindo do vertice ",T.Arc[a].u)
	       	  end
	       	  ub_mtz[T.Arc[a].u] = a
	       end
	   end

	   mtz_cost = SolutionCost(T,ub_mtz)
	   println("Solucao do ciclo obtido pelo MTZ foi de: ",mtz_cost)

	   # --------------------------------------------------
	   if (T.NNodes > maximum_number_of_nodes_to_print)
	      println("Verbose mode: Solution is printed only for graphs up to ",maximum_number_of_nodes_to_print," nodes.")
	   else
		PrintHeader("Solution of the model (variables):")
   	   	for u in 1:T.NNodes
      	       	    println("y_",u," = ",value(y[u]))
   	   	end
   	   	for a in 1:T.NArcs
               	    x_a = value(x[a])
	       	    if (x_a!=0)
	       	       println("x_",T.Arc[a].u,"_",T.Arc[a].v," = ",x_a)
	       	    end
	   	end
	   end
	   
	   # --------------------------------------------------
	   PrintHeader("Solution as a sequence of nodes of the TSP cycle")
	   u = 1
	   for i in 1:T.NNodes
	       print(u,", ")
	       u = T.Arc[ub_mtz[u]].v
	   end
	   println()

	   # --------------------------------------------------
	   PrintHeader("Solution as a sequence of arcs of the TSP cycle (arc_id,u,v)")
	   u = 1
	   for i in 1:T.NNodes
	       print("(",ub_mtz[u],",",T.Arc[ub_mtz[u]].u,",",T.Arc[ub_mtz[u]].v,") ")
	       u = T.Arc[ub_mtz[u]].v
	   end
	   println()

	   # --------------------------------------------------
	   PrintHeader("Cycle containing node 1")
	   u = 1
	   cyclesize=0
	   while (true)
	   	 print(u,", ")
	    	 u = T.Arc[ub_mtz[u]].v
	    	 cyclesize += 1
	    	 if (u==1)
	       	    break
	    	 end
	   end
	   println()
	   println("Node 1 is in a cycle of size ",cyclesize)
	   PrintHeader("")
      end
end
