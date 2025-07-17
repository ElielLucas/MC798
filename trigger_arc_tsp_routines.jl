#
# Este deve ser o unico arquivo que deve ser atualizado pelo aluno
#
# Abaixo voce encontra as rotinas vazias das seguintes funcoes:
#                 TriggerArcTSP_lb_lp(T)
#                 TriggerArcTSP_lb_rlxlag(T)
#                 TriggerArcTSP_lb_colgen(T) - Opcional
#                 TriggerArcTSP_ub_lp(T)
#                 TriggerArcTSP_ub_rlxlag(T)
#                 TriggerArcTSP_ub_colgen(T) - Opcional
#                 TriggerArcTSP_ilp(T)
#
using Dates
using JuMP
using Gurobi
using Printf
using Plots
using Infiltrator
using Graphs
using SimpleWeightedGraphs



function calcular_custo_tatsp(tour_indices::Vector{Int}, arcos::Vector{ArcType}, relacoes::Vector{TriggerType}, num_nos::Int)
    # Reconstrói a ordem dos arcos no ciclo com base na sequência
    tour_arcos = Vector{ArcType}(undef, num_nos)
    visitado = falses(num_nos)
    mapa_saida = Dict{Int, ArcType}()

    for idx in tour_indices
        a = arcos[idx]
        mapa_saida[a.u] = a
    end

    # Reconstruir ciclo ordenado
    atual = 1
    for i in 1:num_nos
        if !haskey(mapa_saida, atual)
            error("Tour incompleto ou inválido. Nó $atual não possui arco de saída.")
        end
        tour_arcos[i] = mapa_saida[atual]
        atual = mapa_saida[atual].v
    end

    # Prepara matriz de gatilhos
    gatilhos = Dict{Tuple{Int, Int}, Float64}()
    for r in relacoes
        gatilhos[(r.trigger_arc_id, r.target_arc_id)] = r.cost
    end

    # Calcula custo com gatilhos
    custo_total = 0.0
    for i in 1:num_nos
        arc_i = tour_indices[i]
        custo_arc_i = arcos[arc_i].cost

        # Verifica se algum arco anterior o ativa
        for j in 1:i-1
            arc_j = tour_indices[j]
            if haskey(gatilhos, (arc_j, arc_i))
                custo_arc_i = gatilhos[(arc_j, arc_i)]
                break
            end
        end

        custo_total += custo_arc_i
    end

    return custo_total
end

# --------------------------------------------------------------
function modelar_tatsp_ilp(n::Int, m::Int, R::Int,
                       arcos::Vector{ArcType},
                       rels::Vector{TriggerType})
    model = Model(Gurobi.Optimizer)
    M = n

    # variáveis
    @variable(model, x[1:m],       Bin)      # seleciona arco
    @variable(model, ya[1:m],      Bin)      # paga custo base
    @variable(model, yr[1:R],      Bin)      # ativa relação
    @variable(model, yhat[1:R],    Bin)      # último gatilho antes do alvo
    @variable(model, 1 <= u[1:n] <= n, Int)  # MTZ

    # (1) objetivo
    @objective(model, Min,
        sum(rels[r].cost * yr[r] for r in 1:R) +
        sum(arcos[a].cost * ya[a] for a in 1:m)
    )

    # (2) escolhe exatamente n arcos
    @constraint(model, sum(x[a] for a in 1:m) == n)

    # (3) subtour‐elimination (MTZ)
    for i in 1:m
        a = arcos[i]
        if a.v != 1
            @constraint(model,
                u[a.u] + 1 <= u[a.v] + M * (1 - x[i])
            )
        end
    end

    # (4) grau de entrada = 1
    for j in 1:n
        @constraint(model,
            sum(x[i] for i in 1:m if arcos[i].v == j) == 1
        )
    end

    # (5) grau de saída = 1
    for i in 1:n
        @constraint(model,
            sum(x[j] for j in 1:m if arcos[j].u == i) == 1
        )
    end

    # (6) vínculo custo base / relação
    for i in 1:m
        Ra = [r_idx for (r_idx, r) in enumerate(rels) if r.target_arc_id == i]
        @constraint(model,
            ya[i] + sum(yr[r] for r in Ra) == x[i]
        )
    end

    # (7)–(10) ativação de relações
    for (r_idx, r) in enumerate(rels)
        t = arcos[r.trigger_arc_id]
        a = arcos[r.target_arc_id]
        @constraint(model, yr[r_idx] <= x[r.trigger_arc_id])                                
        @constraint(model,
            u[t.u] + 1 <= u[a.u] + M * (1 - yr[r_idx])
        )
        @constraint(model,
            u[a.u] + 1 <= u[t.u] + M * (1 - yhat[r_idx])
        )
        @constraint(model,
            x[r.trigger_arc_id]
            <= (1 - x[r.target_arc_id]) + (1 - ya[r.target_arc_id]) + yhat[r_idx]
        )
    end

    # (11) último gatilho antes do alvo
    for i in 1:R, j in 1:R
        r1, r2 = rels[i], rels[j]
        if r1.target_arc_id == r2.target_arc_id && i != j
            t1 = arcos[r1.trigger_arc_id]
            t2 = arcos[r2.trigger_arc_id]
            @constraint(model,
                u[t1.u] >= u[t2.u]
                            - M * (2 - yr[i] - (x[r2.trigger_arc_id] - yhat[j]))
            )
        end
    end

    return model
end

function resolver_subproblema_deg(λ_out, λ_in, n, m, R, arcos::Vector{ArcType}, rels::Vector{TriggerType})
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    M = n + 1

    @variable(model, 0 <= x[1:m] <= 1)
    @variable(model, 0 <= ya[1:m] <= 1)
    @variable(model, 0 <= yr[1:R] <= 1)
    @variable(model, 0 <= yhat[1:R] <= 1)
    @variable(model, 1 <= u[1:n] <= n)

    # (1) Função objetivo
    @objective(model, Min,
        sum(rels[j].cost * yr[j] for j in 1:R) +
        sum(arcos[i].cost * ya[i] for i in 1:m) +
        sum(λ_out[arcos[i].u] * x[i] for i in 1:m) +
        sum(λ_in[arcos[i].v] * x[i] for i in 1:m)
    )

    # (2) Seleção de |N| arcos
    @constraint(model, sum(x[a] for a in 1:m) == n)

    # (3) Subtour-elimination (MTZ)
    for i in 1:m
        a = arcos[i]
        if a.v != 1
            @constraint(model, u[a.u] + 1 <= u[a.v] + M * (1 - x[i]))
        end
    end

    # (6) Custo base / relação
    for i in 1:m
        Ra = [r_idx for (r_idx, rel) in enumerate(rels) if rel.target_arc_id == i]
        @constraint(model, ya[i] + sum(yr[r] for r in Ra) == x[i])
    end

    # (7)-(10) Ativação de relações
    for (r_idx, r) in enumerate(rels)
        t_idx = r.trigger_arc_id
        a_idx = r.target_arc_id
        t = arcos[t_idx]
        a = arcos[a_idx]

        @constraint(model, yr[r_idx] <= x[t_idx])                                                   # (7)
        @constraint(model, u[t.u] + 1 <= u[a.u] + M*(1 - yr[r_idx]))                                # (8)
        @constraint(model, u[a.u] + 1 <= u[t.u] + M*(1 - yhat[r_idx]))                              # (9)
        @constraint(model, x[t_idx] <= (1 - x[a_idx]) + (1 - ya[a_idx]) + yhat[r_idx])              # (10)
    end

    # (11) Último trigger antes de alvo
    for r1_idx in 1:R, r2_idx in 1:R
        r1 = rels[r1_idx]
        r2 = rels[r2_idx]
        if r1.target_arc_id == r2.target_arc_id && r1_idx != r2_idx
            t1 = arcos[r1.trigger_arc_id]
            t2 = arcos[r2.trigger_arc_id]
            @constraint(model,
                u[t1.u] >= u[t2.u] - M*(2 - yr[r1_idx] - (x[r2.trigger_arc_id] - yhat[r2_idx]))
            )
        end
    end

    optimize!(model)
    return value.(x), objective_value(model)
end


function dual_lagrange_deg(n, m, R, arcos, rels, num_iter, ub_ref, max_time)
    global λ_out = zeros(n)
    global λ_in = zeros(n)
    best_lb = -Inf
    hist = Float64[]
    α = 1.0
    best_x = nothing
    start_time = time()

    for k in 1:num_iter
        if time() - start_time > max_time
            println(@sprintf("Tempo máximo de %.1fs excedido na relaxação Lagrangeana. Interrompendo.", max_time))
            break
        end

        x_sol, Z = resolver_subproblema_deg(λ_out, λ_in, n, m, R, arcos, rels)
        if Z > best_lb
            best_lb = Z
            best_x = copy(x_sol)
        end
        push!(hist, best_lb)

        g_out = zeros(n)
        g_in = zeros(n)
        for i in 1:m
            a = arcos[i]
            g_out[a.u] += x_sol[i]
            g_in[a.v] += x_sol[i]
        end
        g_out .-= 1
        g_in .-= 1

        sq_norm = sum(g_out.^2) + sum(g_in.^2)
        if sq_norm < 1e-8
            break
        end

        t = α * (ub_ref - Z) / sq_norm
        t = clamp(t, -100.0, 100.0)

        λ_out .+= t .* g_out
        λ_in .+= t .* g_in

        α *= 0.95

        println(@sprintf("It %3d: LB=%.2f (best=%.2f)  step=%.2f", k, Z, best_lb, t))
    end

    return best_lb, hist, best_x
end


function solve_lp_relaxation(n::Int, m::Int, R::Int, arcos::Vector{ArcType}, rels::Vector{TriggerType},
                             fixar_uns::Vector{Int} = Int[], fixar_zeros::Vector{Int} = Int[])
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    M = n

    @variable(model, 0 <= x[1:m] <= 1)      # x_a ∀ a ∈ A
    @variable(model, 0 <= ya[1:m] <= 1)     # y_a ∀ a ∈ A
    @variable(model, 0 <= yr[1:R] <= 1)     # y_r ∀ r ∈ R
    @variable(model, 0 <= yhat[1:R] <= 1)   # ŷ_r ∀ r ∈ R
    @variable(model, 1 <= u[1:n] <= n, Int) # u_i ∈ {1…n}

    for i in fixar_uns
        @constraint(model, x[i] == 1)
    end

    for i in fixar_zeros
        @constraint(model, x[i] == 0)
    end

    # (1) Função-objetivo
    @objective(model, Min,
        sum(rels[r].cost * yr[r] for r in 1:R) +
        sum(arcos[a].cost * ya[a] for a in 1:m)
    )

    # (2) Seleção de |N| arcos
    @constraint(model, sum(x[a] for a in 1:m) == n)

    # (3) Subtour-elimination (MTZ)
    for a in 1:m
        arc = arcos[a]
        if arc.v != 1
            @constraint(model, u[arc.u] + 1 <= u[arc.v] + M * (1 - x[a]))
        end
    end

    # (4) Grau de entrada
    for j in 1:n
        @constraint(model,
            sum(x[a] for a in 1:m if arcos[a].v == j) == 1
        )
    end

    # (5) Grau de saída
    for i in 1:n
        @constraint(model,
            sum(x[a] for a in 1:m if arcos[a].u == i) == 1
        )
    end

    # (6) Custo base / relação
    for a in 1:m
        Ra = [r_idx for (r_idx, r) in enumerate(rels) if r.target_arc_id == a]
        @constraint(model, ya[a] + sum(yr[r] for r in Ra) == x[a])
    end

    # (7)-(10) Ativação de relações
    for (r_idx, r) in enumerate(rels)
        t_idx = r.trigger_arc_id
        a_idx = r.target_arc_id
        t = arcos[t_idx]
        a = arcos[a_idx]

        @constraint(model, yr[r_idx] <= x[t_idx])                                                    # (7)
        @constraint(model, u[t.u] + 1 <= u[a.u] + M * (1 - yr[r_idx]))                               # (8)
        @constraint(model, u[a.u] + 1 <= u[t.u] + M * (1 - yhat[r_idx]))                             # (9)
        @constraint(model, x[t_idx] <= (1 - x[a_idx]) + (1 - ya[a_idx]) + yhat[r_idx])               # (10)
    end

    # (11) Último trigger antes de alvo
    for r1_idx in 1:R, r2_idx in 1:R
        r1 = rels[r1_idx]
        r2 = rels[r2_idx]
        if r1.target_arc_id == r2.target_arc_id && r1_idx != r2_idx
            t1 = arcos[r1.trigger_arc_id]
            t2 = arcos[r2.trigger_arc_id]
            @constraint(model,
                u[t1.u] >= u[t2.u] - M * (2 - yr[r1_idx] - (x[r2.trigger_arc_id] - yhat[r2_idx]))
            )
        end
    end

    t0 = time()
    optimize!(model)
    lp_time = time() - t0
    lb = objective_value(model)
    println(@sprintf("LP solved in %.3fs, LB=%.2f", lp_time, lb))

    return value.(x), lp_time, lb
end

function three_opt!(T::TriggerArcTSP, tour::Vector{Int})
    best_cost = calcular_custo_tatsp(tour, T.Arc, T.Trigger, T.NNodes)

    order = [ T.Arc[i].u for i in tour ]
    push!(order, T.Arc[last(tour)].v)
    nn = length(order)

    for i in 2:nn-3
      for j in i+1:nn-2
        for k in j+1:nn-1
          seg1 = order[1:i]
          seg2 = order[i+1:j]
          seg3 = order[j+1:k]
          seg4 = order[k+1:end]

          orders = [
            vcat(seg1, reverse(seg2), seg3, seg4),
            vcat(seg1, seg2, reverse(seg3), seg4),
            vcat(seg1, reverse(seg2), reverse(seg3), seg4),
            vcat(seg1, seg3, seg2, seg4),
            vcat(seg1, seg3, reverse(seg2), seg4),
            vcat(seg1, reverse(seg3), seg2, seg4),
            vcat(seg1, reverse(seg3), reverse(seg2), seg4)
          ]

          for new_order in orders
            new_tour = Int[]
            valid = true
            for p in 1:(length(new_order)-1)
              o, d = new_order[p], new_order[p+1]
              f = findfirst(a -> a.u == o && a.v == d, T.Arc)
              if f === nothing
                valid = false
                break
              end
              push!(new_tour, f)
            end

            if valid && length(new_tour) == length(tour)
              c = calcular_custo_tatsp(new_tour, T.Arc, T.Trigger, T.NNodes)
              if c < best_cost
                tour[:] = new_tour
                return true
              end
            end
          end
        end
      end
    end

    return false
end

function run_rounding(T::TriggerArcTSP, x_frac::Vector{Float64};
                     threshold=0.5, k=3, num_iter=100)

    n = T.NNodes
    best_cost = Inf
    times = Float64[]
    costs = Float64[]
    elapsed = 0.0
    t_start = time()

    for iter in 1:num_iter
        # critério de parada por tempo
        if time() - t_start > T.maxtime_ub_lp
            println(@sprintf("Rounding: tempo máximo de %ds excedido, saindo no it %d.", 
                             T.maxtime_ub_lp, iter))
            break
        end

        t1 = time()

        tour, cost0 = greedy_guided_construction(T, x_frac)
        
        if !isempty(tour)
            three_opt!(T, tour)
            cost1 = calcular_custo_tatsp(tour, T.Arc, T.Trigger, T.NNodes)
        else
            cost1 = cost0
        end
        
        best_cost = min(best_cost, cost1)
        push!(times, elapsed)
        push!(costs, best_cost)
        elapsed += time() - t1
        
        println(@sprintf("Iter %3d: best=%.2f @ t=%.3f", iter, best_cost, elapsed))
    end

    return best_cost, times, costs
end

function tour_to_node_sequence(T::TriggerArcTSP, tour_arcs::Vector{Int})
    node_seq = [1]               
    for a in tour_arcs
        push!(node_seq, T.Arc[a].v)
    end
    return node_seq
end

function node_sequence_to_arcs(T::TriggerArcTSP, node_seq::Vector{Int})
    tour = Int[]
    for i in 1:length(node_seq)-1
        u, v = node_seq[i], node_seq[i+1]
        f = findfirst(a -> a.u == u && a.v == v, T.Arc)
        if f === nothing
            error("Tour incompleto: não existe arco $u → $v")
        end
        push!(tour, f)
    end
    return tour
end

function perturbar_seq_nos_conectada!(T::TriggerArcTSP, seq::Vector{Int})
    valid_arc = Dict((a.u, a.v) => true for a in T.Arc)

    for _ in 1:50
        i, j = sort(rand(2:length(seq)-1, 2))
        if j - i < 1
            continue
        end

        cand = copy(seq)
        cand[i:j] = reverse(cand[i:j])

        valid = all(haskey(valid_arc, (cand[k], cand[k+1])) for k in 1:length(cand)-1)

        if valid
            seq[i:j] = reverse(seq[i:j])
            return true
        end
    end

    return false
end

function perturbar_seq_nos_conectada_forte!(T::TriggerArcTSP, seq::Vector{Int})
    valid_arc = Dict((a.u, a.v) => true for a in T.Arc)

    for _ in 1:50
        i, j = sort(rand(2:length(seq)-1, 2))
        k, l = sort(rand(2:length(seq)-1, 2))
        if j - i < 1 || l - k < 1 || i == k
            continue
        end

        seg1 = seq[i:j]
        seg2 = seq[k:l]

        cand = copy(seq)
        cand[i:j] = seg2[1:min(end, j-i+1)]
        cand[k:l] = seg1[1:min(end, l-k+1)]

        valid = all(haskey(valid_arc, (cand[q], cand[q+1])) for q in 1:length(cand)-1)

        if valid
            seq .= cand
            return true
        end
    end

    return false
end

function construir_tour_insercao_xfrac(x_frac::Vector{Float64},
                                        arcos::Vector{ArcType},
                                        n::Int)
    arcs_from_node = [Int[] for _ in 1:n]
    for (i, a) in enumerate(arcos)
        push!(arcs_from_node[a.u], i)
    end

    visited = falses(n)
    current = 1
    visited[current] = true

    tour = Int[]
    for step in 1:n
        best_idx = 0
        best_score = Inf
        is_last = (length(tour) == n - 1)

        for i in arcs_from_node[current]
            a = arcos[i]
            dest = a.v

            if (!visited[dest] && !is_last) || (is_last && dest == 1)
                score = a.cost - 1000.0 * x_frac[i]
                if score < best_score
                    best_score = score
                    best_idx = i
                end
            end
        end

        if best_idx == 0
            println("Inserção: sem arco viável em nó $current → tour incompleto.")
            break
        end

        push!(tour, best_idx)
        current = arcos[best_idx].v
        visited[current] = true
    end

    if length(tour) == n && current != 1
        println("Inserção: tour tem $n arcos mas não retorna ao depósito.")
    elseif length(tour) != n
        println("Inserção: tour gerou $(length(tour)) arcos, esperado $n.")
    end

    return tour
end


function bfs_to_nearest_unvisited(T::TriggerArcTSP,
                                  start::Int,
                                  visited::Set{Int})
    queue = [(start, Int[])]
    seen  = Set([start])
    while !isempty(queue)
        u, path = popfirst!(queue)
        for (i,a) in enumerate(T.Arc)
            if a.u == u && !(a.v in seen)
                new_path = copy(path)
                push!(new_path, i)
                if !(a.v in visited)
                    return new_path, a.v
                end
                push!(queue, (a.v, new_path))
                push!(seen, a.v)
            end
        end
    end
    error("bfs_to_nearest_unvisited: não há caminho para nenhum não-visitado")
end

function greedy_complete_tour_xfrac_bfs(
        T::TriggerArcTSP,
        xfrac::Vector{Float64};
        α::Float64=1000.0
    )
    n       = T.NNodes
    visited = Set([1])
    tour    = Int[]
    current = 1

    while length(visited) < n
        cands = [(a.v, i, a.cost, xfrac[i]) for (i,a) in enumerate(T.Arc)
                 if a.u == current && !(a.v in visited)]
        if !isempty(cands)
            scores = [c[3] - α*c[4] for c in cands]
            idx    = argmin(scores)
            v, arc_idx, _, _ = cands[idx]
            push!(tour, arc_idx)
            current = v
            push!(visited, v)
        else
            path_arcs, new_node = bfs_to_nearest_unvisited(T, current, visited)
            for aidx in path_arcs
                push!(tour, aidx)
                current = T.Arc[aidx].v
                push!(visited, current)
            end
        end
    end

    back = findfirst(a -> a.u == current && a.v == 1, T.Arc)
    if back !== nothing
        push!(tour, back)
    else
        path_arcs, _ = bfs_to_nearest_unvisited(T, current, Set([1]))
        append!(tour, path_arcs)
    end

    return tour
end


function greedy_guided_construction(
        T::TriggerArcTSP,
        best_x::Vector{Float64};
        α::Float64 = 1000.0
    )
    tour = greedy_complete_tour_xfrac_bfs(T, best_x; α=α)
    cost = calcular_custo_tatsp(tour, T.Arc, T.Trigger, T.NNodes)

    return tour, cost
end


function ils_3opt!(
    T::TriggerArcTSP,
    xfrac::Vector{Float64};
    max_iter=1000,
    max_no_improve=100,
    α_accept=1.05,
    restart_every=200
)
    n = T.NNodes
    t_start = time()

    tour, _ = greedy_guided_construction(T, xfrac)
    while three_opt!(T, tour); end

    best_tour = copy(tour)
    best_seq  = tour_to_node_sequence(T, best_tour)
    best_cost = calcular_custo_tatsp(best_tour, T.Arc, T.Trigger, n)

    current_tour = copy(best_tour)
    current_seq  = copy(best_seq)
    current_cost = best_cost

    no_improve = 0
    hist = [best_cost]

    println(@sprintf("ILS-3OPT: custo inicial = %.2f", best_cost))
    for it in 1:max_iter
        if time() - t_start > T.maxtime_ub_rlxlag
            println("Tempo máximo excedido. Saindo.")
            break
        end

        if restart_every > 0 && it % restart_every == 0
            println(">> Reinício aleatório")
            current_tour, _ = greedy_guided_construction(T, rand(length(xfrac)))
            while three_opt!(T, current_tour); end
            current_seq  = tour_to_node_sequence(T, current_tour)
            current_cost = calcular_custo_tatsp(current_tour, T.Arc, T.Trigger, n)
        end

        seq = copy(current_seq)
        perturbar_seq_nos_conectada!(T, seq)

        tour = node_sequence_to_arcs(T, seq)
        while three_opt!(T, tour); end

        cost = calcular_custo_tatsp(tour, T.Arc, T.Trigger, n)

        if cost < best_cost
            best_cost = cost
            best_tour = copy(tour)
            best_seq  = copy(seq)
            println("Melhoria absoluta: $cost")
            no_improve = 0
        elseif cost < α_accept * best_cost
            current_tour = copy(tour)
            current_seq  = copy(seq)
            current_cost = cost
            println("Aceita pior solução com custo $cost")
            no_improve += 1
        else
            no_improve += 1
        end

        push!(hist, best_cost)
        println(@sprintf("ILS-3OPT: It %3d: Custo = %.2f | Melhor = %.2f", it, cost, best_cost))

        if no_improve >= max_no_improve
            println("Sem melhora por $max_no_improve iterações.")
            break
        end
    end

    return best_tour, best_cost, hist
end


function find_violated_subtour(T::TriggerArcTSP, x_frac::Vector{Float64}; tol=1e-6)
    g = SimpleDiGraph(T.NNodes)
    for (i, a) in enumerate(T.Arc)
        if x_frac[i] > tol
            add_edge!(g, a.u, a.v)
        end
    end

    for S in strongly_connected_components(g)
        if 1 < length(S) < T.NNodes
            return S
        end
    end
    return nothing
end

function plot_and_save(x, y; lb=nothing, title="", xlabel="", ylabel="", filename="", legend_label="", hline_label="")
    p = plot(x, y,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        label=legend_label == "" ? nothing : legend_label,
        marker=:circle,
        legend=legend_label != "")

    if lb !== nothing
        hline!([lb], linestyle=:dash, label=hline_label)
    end

    savefig(p, filename)
    println("Gráfico salvo: $filename")
end
# --------------------------------------------------------------

function TriggerArcTSP_lb_lp(T::TriggerArcTSP)
	println("---- Iniciando Relaxamento LP ----")
    Random.seed!(T.seednumber)
    StartingTime = time()
    
	n, m, R = T.NNodes, T.NArcs, T.NTriggers

    x_frac, lp_time, lb = solve_lp_relaxation(n, m, R, T.Arc, T.Trigger)
    T.lb_lp = lb
    
    T.time_lb_lp = ceil(Int, time() - StartingTime)
	global x_frac_lb_lp = x_frac
end
# --------------------------------------------------------------

function TriggerArcTSP_lb_rlxlag(T::TriggerArcTSP)
	println("---- Iniciando Relaxação Lagrangeana ----")
    Random.seed!(T.seednumber)
    StartingTime = time()

    n = T.NNodes
    m = T.NArcs
    R = T.NTriggers

    arcos = T.Arc
    rels = T.Trigger

	best_cost, _, _ = run_rounding(T, x_frac_lb_lp; threshold=0.5, k=3, num_iter=100)
    global ub_ref = best_cost

	lb, hist_lb, best_x = dual_lagrange_deg(n, m, R, arcos, rels, 1000, ub_ref, T.maxtime_lb_rlxlag)

    T.lb_rlxlag = lb
    
	T.time_lb_rlxlag = ceil(Int, time() - StartingTime)

	plot_and_save(1:length(hist_lb), hist_lb;
	    lb=maximum(hist_lb),
	    title="Relaxação Lagrangeana",
	    xlabel="Iteração", ylabel="Limitante Inferior",
	    filename="lb_rlxlag_$(T.inputfilename[1:end-4]).png",
	    legend_label="Melhor LB", hline_label="Melhor LB Final")
	
    global best_x_rlxlag = best_x
	# TriggerArcTSP_lb_fix(T)
end
# --------------------------------------------------------------
function TriggerArcTSP_lb_colgen(T::TriggerArcTSP)
    println("---- Iniciando Lower Bound via Column Generation (tours completos) ----")
    Random.seed!(T.seednumber)
    StartingTime = time()

    n, m, R = T.NNodes, T.NArcs, T.NTriggers

    # 1) Garante solução Lagrangeana disponível
    if !@isdefined(best_x_rlxlag)
        TriggerArcTSP_lb_rlxlag(T)
    end

    # 2) Gera várias colunas iniciais (tours heurísticos distintos)
    columns_arcs  = Vector{Vector{Int}}()
    columns_nodes = Vector{BitVector}()
    col_costs     = Float64[]

    # 2.1) uma via greedy_guided_construction
    arc_seq, cost = greedy_guided_construction(T, best_x_rlxlag)
    if !isempty(arc_seq)
        push!(columns_arcs, arc_seq)
        bv = falses(n)
        for v in tour_to_node_sequence(T, arc_seq)[1:end-1]
            bv[v] = true
        end
        push!(columns_nodes, bv)
        push!(col_costs, cost)
    end

    # 2.2) mais algumas via greedy_guided_construction com ruído
    for rep in 1:5
        xf = [ x + 0.05*randn() for x in best_x_rlxlag ]
        seq, _ = greedy_guided_construction(T, xf)
        if length(seq)==n && !(any(a->a==seq, columns_arcs))
            c = calcular_custo_tatsp(seq, T.Arc, T.Trigger, T.NNodes)
            push!(columns_arcs, seq)
            bv = falses(n)
            for v in tour_to_node_sequence(T, seq)[1:end-1]
                bv[v] = true
            end
            push!(columns_nodes, bv)
            push!(col_costs, c)
        end
    end

    # 3) loop de column-generation
    for iter in 1:30
        # 3a) resolve mestre restrito θ ≥ 0
        master = Model(Gurobi.Optimizer); set_silent(master)
        k = length(columns_arcs)
        @variable(master, θ[1:k] >= 0)
        cons = Vector{ConstraintRef}(undef, n)
        for j in 1:n
            cons[j] = @constraint(master,
                sum(columns_nodes[c][j] ? θ[c] : 0 for c in 1:k) == 1
            )
        end
        @objective(master, Min, sum(col_costs[c]*θ[c] for c in 1:k))
        optimize!(master)
        if termination_status(master) != MOI.OPTIMAL
            error("LB_colgen: mestre restrito não ótimo")
        end
        lb = objective_value(master)
        println(@sprintf("  it %2d → LB mestre = %.4f", iter, lb))

        # 3b) extrai duais π_j
        π = [ dual(cons[j]) for j in 1:n ]

        # 3c) pricing exato: resolve TSP com custos reduzidos
        pricing = modelar_tatsp_ilp(n, m, R, T.Arc, T.Trigger)
        @objective(pricing, Min,
            sum((T.Arc[a].cost + π[T.Arc[a].u] + π[T.Arc[a].v]) * pricing[:x][a]
                for a in 1:m)
        )
        optimize!(pricing)
        if termination_status(pricing) != MOI.OPTIMAL
            error("LB_colgen: pricing MIP não ótimo")
        end

        # 3d) extrai tour binário do pricing
        xbin = [ round(value(pricing[:x][a])) for a in 1:m ]
        new_arc_seq = findall(v->v==1, xbin)
        if length(new_arc_seq) != n
            println("    → pricing não gerou tour hamiltoniano; interrompendo.")
            break
        end

        # 3e) calcula reduced cost e decide parar/adicionar
        cost_new = calcular_custo_tatsp(new_arc_seq, T.Arc, T.Trigger, T.NNodes)
        rc = objective_value(pricing) - sum(π)
        println(@sprintf("    → rc = %.4f  (custo=%.2f)", rc, cost_new))
        if rc >= -1e-6
            println("    → convergiu; LB_final = $lb")
            T.lb_colgen      = lb
            T.time_lb_colgen = ceil(time() - StartingTime)
            global paths_colgen_arcs, paths_colgen_nodes, col_costs_colgen
            paths_colgen_arcs  = columns_arcs
            paths_colgen_nodes = columns_nodes
            col_costs_colgen   = col_costs
            return
        end

        # 3f) adiciona nova coluna
        bv = falses(n)
        for v in tour_to_node_sequence(T, new_arc_seq)[1:end-1]
            bv[v] = true
        end
        # evita duplicata exata
        if !any(bv == bn for bn in columns_nodes)
            push!(columns_arcs, new_arc_seq)
            push!(columns_nodes, bv)
            push!(col_costs, cost_new)
        else
            println("    → coluna repetida, parando.")
            break
        end
    end

    # não convergiu em 30 iterações
    println("LB_colgen: max it atingido → LB≈$(col_costs[end])")
    T.lb_colgen      = col_costs[end]
    T.time_lb_colgen = ceil(time() - StartingTime)
    global paths_colgen_arcs, paths_colgen_nodes, col_costs_colgen
    paths_colgen_arcs  = columns_arcs
    paths_colgen_nodes = columns_nodes
    col_costs_colgen   = col_costs
end

# --------------------------------------------------------------
function TriggerArcTSP_ub_lp(T::TriggerArcTSP)
	println("---- Iniciando Heurística de Arredondamento ----")
    Random.seed!(T.seednumber)
	StartingTime = time()

	best_cost, times, costs = run_rounding(T, x_frac_lb_lp; threshold=0.5, k=3, num_iter=1000)

	T.time_ub_lp = ceil(time() - StartingTime)

	plot_and_save(times, costs;
        lb=T.lb_lp,
        title="Convergência da Heurística de Arredondamento",
        xlabel="Tempo (s)", ylabel="Melhor Custo",
        filename="ub_lp_$(T.inputfilename[1:end-4]).png",
        legend_label="", hline_label="LP Lower Bound")

	T.ub_lp = best_cost
end

# --------------------------------------------------------------
function TriggerArcTSP_ub_rlxlag(T::TriggerArcTSP)
	println("---- Iniciando ILS completo ----")
    Random.seed!(T.seednumber)
    StartingTime = time()

    best_x = best_x_rlxlag
    tour, cost, hist = ils_3opt!(T, best_x;
                                  max_iter=1000,
                                  max_no_improve=200,
                                  α_accept=1.05,
                                  restart_every=200)

    T.ub_rlxlag_arcs = tour
    T.ub_rlxlag      = cost
    T.time_ub_rlxlag = ceil(Int, time() - StartingTime)

	plot_and_save(1:length(hist), hist;
        lb=minimum(hist),
        title="ILS com 3-Opt",
        xlabel="Iteração", ylabel="Custo",
        filename="ub_rlxlag_$(T.inputfilename[1:end-4]).png",
        legend_label="Melhor Custo", hline_label="Melhor Custo Final")
end
# --------------------------------------------------------------
function TriggerArcTSP_ub_colgen(T::TriggerArcTSP)
    println("---- Iniciando Upper Bound via Column Generation ----")
    Random.seed!(T.seednumber)
    StartingTime = time()

    # garante LB já rodou
    if !@isdefined(paths_colgen_arcs)
        TriggerArcTSP_lb_colgen(T)
    end

    columns_arcs = paths_colgen_arcs
    columns_nodes = paths_colgen_nodes
    costs = col_costs_colgen
    n = T.NNodes
    k = length(columns_arcs)

    master = Model(Gurobi.Optimizer); set_silent(master)
    @variable(master, θ[1:k], Bin)
    for j in 1:n
        @constraint(master,
            sum(columns_nodes[c][j] ? θ[c] : 0 for c in 1:k) == 1
        )
    end
    @objective(master, Min, sum(costs[c]*θ[c] for c in 1:k))
    optimize!(master)
    if termination_status(master) != MOI.OPTIMAL
        error("UB_colgen: mestre MIP não ótimo")
    end

    vals = value.(θ)
    idx  = findfirst(v->isapprox(v,1; atol=1e-6), vals)
    T.ub_colgen      = costs[idx]
    T.time_ub_colgen = ceil(time() - StartingTime)
    println(@sprintf("    → UB_final = %.2f", T.ub_colgen))
end

# --------------------------------------------------------------

function TriggerArcTSP_ilp(T::TriggerArcTSP)
    println("---- Iniciando Algoritmo Exato - PLI com Gurobi ----")
    start_time = time()
    Random.seed!(T.seednumber)

    n, m, R = T.NNodes, T.NArcs, T.NTriggers
    arcos, rels = T.Arc, T.Trigger

    model = modelar_tatsp_ilp(n, m, R, arcos, rels)
    set_optimizer_attribute(model, "TimeLimit", T.maxtime_ilp)

    optimize!(model)

    if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.TIME_LIMIT
        x_val = value.(model[:x])
        tour_indices = findall(v -> isapprox(v, 1.0; atol=1e-4), x_val)

        if length(tour_indices) != n
            println("!!! Tour inválido: esperados $n arcos, obtidos $(length(tour_indices)).")
        end

        custo_real = calcular_custo_tatsp(tour_indices, arcos, rels, n)

        T.ub_ilp = custo_real
        T.lb_ilp = objective_value(model)
        T.ub_ilp_arcs = tour_indices
        T.nn_ilp = MOI.get(model, MOI.NodeCount())
    else
        println("Erro: Solver não encontrou solução viável.")
    end

    T.time_ilp = ceil(Int, time() - start_time)
    println("---- PLI com Gurobi Finalizado ----")
    println("Custo ótimo (real): ", T.ub_ilp)
    println("LB solver (função objetivo): ", T.lb_ilp)
    println("Arcos no tour: ", T.ub_ilp_arcs)
end
