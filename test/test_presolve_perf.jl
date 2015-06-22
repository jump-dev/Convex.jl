using Convex
import Distributions
using FactCheck
using SCS

facts("Presolve performance") do
  context("svm") do
    function gen_data(n)
      pos = rand(Distributions.MvNormal([1.0,2.0],1.0),n)
      neg = rand(Distributions.MvNormal([-1.0,1.0],1.0),n)
      return pos,neg
    end

    const N = 2
    const C = 10

    function svm(pos, neg, use_presolve)
      w = Variable(N)
      b = Variable()
      ξpos = Variable(size(pos,2), Positive())
      ξneg = Variable(size(neg,2), Positive())

      problem = minimize(sum_squares(w) + C*sum(ξpos) + C*sum(ξneg))
      for j in 1:size(pos,2)
        push!(problem.constraints, dot(w,pos[:,j]) - b >= 1 - ξpos[j])
        #problem.constraints += dot(w,pos[:,j]) - b >= 1 - ξpos[j]
      end
      for j in 1:size(neg,2)
        push!(problem.constraints, -1*(dot(w,neg[:,j]) - b) >= 1-ξneg[j])
        #problem.constraints += -1*(dot(w,neg[:,j]) - b) >= 1-ξneg[j]
      end
      solve!(problem, SCSSolver(verbose=0), use_presolve=use_presolve)
      return evaluate(w), evaluate(b)
    end

    # initial compilation
    pos,neg = gen_data(10)
    svm(pos, neg, true)
    # svm(pos, neg, false)

    pos,neg = gen_data(2000)
    @time w, b = svm(pos, neg, true)
    @show w,b
    @time w, b = svm(pos, neg, false)
    @show w,b
  end

  context("variables to constants") do
    # warm up
    x = Variable()
    p = minimize(sum(x), [x == 1])
    solve!(p, SCSSolver(verbose=0))

    x1 = Variable(1000)
    p1 = minimize(sum(x1), [x1 == [1:1000]])
    @time solve!(p1, SCSSolver(verbose=0), use_presolve=false)

    x2 = Variable(1000)
    p2 = minimize(sum(x2), [x2 == [1:1000]])
    @time solve!(p2, SCSSolver(verbose=0), use_presolve=true)
  end
end







