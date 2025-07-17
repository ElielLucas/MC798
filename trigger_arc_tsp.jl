#--------------------------------------------------------
#        Nao mude o conteudo deste arquivo
#--------------------------------------------------------
using ArgParse
using Base
using Random

#--------------------------------------------------------
# Mude apenas o conteudo do arquivo trigger_arc_tsp_routines.jl

include("trigger_arc_tsp_include.jl")	# Nao mude este arquivo

include("trigger_arc_tsp_routines.jl")	# Coloque suas rotinas neste arquivo

include("exemplo_usando_pli.jl")

function main()
	 
	 (inputfilename,	# input file as given in the internet
	 seednumber,		# Randomized steps become deterministic for fixed seed
	 maxtime_lb_lp,		# Max. time to compute LB via LP (may use cutting planes)
	 maxtime_lb_rlxlag,	# Max. time to compute LB by Lagrangean Relaxation
	 maxtime_lb_colgen,	# Max. time to compute LB via LP column generation
	 maxtime_ub_lp,		# Max. time to compute UB via LP rounding
	 maxtime_ub_rlxlag,	# Max. time to compute UB via Lagrangean Relaxation
	 maxtime_ub_colgen,  	# Max. time to compute UB via LP column generation
	 maxtime_ilp,		# Max. time to compute LB and UB via exact Branch and Cut
	 ra,logfilename) = getparameters()

	 T=TriggerArcTSP(inputfilename,seednumber,
			 maxtime_lb_lp,
			 maxtime_lb_rlxlag, 
	 		 maxtime_lb_colgen,
			 maxtime_ub_lp,
			 maxtime_ub_rlxlag,
			 maxtime_ub_colgen,
			 maxtime_ilp,
			 ra,logfilename)

	 # Rotinas desenvolvidas pelos alunos
	 TriggerArcTSP_lb_lp(T)
	 WriteLogFile(T,"lb_lp")
	 
	 TriggerArcTSP_lb_rlxlag(T)
	 WriteLogFile(T,"lb_rlxlag")
	 
	 TriggerArcTSP_ub_lp(T)
	 WriteLogFile(T,"ub_lp")

	 TriggerArcTSP_ub_rlxlag(T)
	 WriteLogFile(T,"ub_rlxlag")

	#  TriggerArcTSP_lb_colgen(T)
	#  WriteLogFile(T,"lb_colgen")

	#  TriggerArcTSP_ub_colgen(T)
	#  WriteLogFile(T,"ub_colgen")
	 
	#  TriggerArcTSP_ilp(T)
	#  WriteLogFile(T,"ilp")

	#  Exemplo_PLI(T,1000000) # The second parameter is the maximum time in seconds
end

# --------------------------------------------------------------------
main()
