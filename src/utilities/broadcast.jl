# This file adds a new syntax to call elementwise operations on Convex expressions
# because overloading broadcast no longer works in Julia v0.6

if VERSION>=v"0.6.0-pre"
  import TakingBroadcastSeriously: unfuse, broadcast_
  export broadcast_

  unfuse(AbstractExpr)

  broadcast_(::typeof(*), x::Constant, y::AbstractExpr) = dot_multiply(x, y)
  broadcast_(::typeof(*), y::AbstractExpr, x::Constant) = dot_multiply(x, y)
  broadcast_(::typeof(*), x::Value, y::AbstractExpr) = dot_multiply(Constant(x), y)
  broadcast_(::typeof(*), x::AbstractExpr, y::Value) = dot_multiply(Constant(y), x)
  broadcast_(::typeof(/), x::AbstractExpr, y::Value) = dot_multiply(Constant(1./y), x)
  broadcast_(::typeof(/), x::Value, y::AbstractExpr) = QolElemAtom(Constant(sqrt(x)), y)
  broadcast_(::typeof(^), x::AbstractExpr, k::Int) = k==2 ? QolElemAtom(x, Constant(ones(size(x)))) : error("raising variables to powers other than 2 is not implemented")
end

if VERSION<v"0.6.0-pre"
  .*(x::Value, y::AbstractExpr) = broadcast(*, x, y)
  .*(x::AbstractExpr, y::Value) = broadcast(*, x, y)
  .*(x::AbstractExpr, y::AbstractExpr) = broadcast(*, x, y)
  ./(x::Value, y::AbstractExpr) = broadcast(/, x, y)
  ./(x::AbstractExpr, y::Value) = broadcast(/, x, y)
  ./(x::AbstractExpr, y::AbstractExpr) = broadcast(/, x, y)
  .^(x::AbstractExpr, y::Value) = broadcast(^, x, y)

  broadcast(::typeof(*), x::Constant, y::AbstractExpr) = dot_multiply(x, y)
  broadcast(::typeof(*), y::AbstractExpr, x::Constant) = dot_multiply(x, y)
  broadcast(::typeof(*), x::Value, y::AbstractExpr) = dot_multiply(Constant(x), y)
  broadcast(::typeof(*), x::AbstractExpr, y::Value) = dot_multiply(Constant(y), x)
  broadcast(::typeof(/), x::AbstractExpr, y::Value) = dot_multiply(Constant(1./y), x)
  broadcast(::typeof(/), x::Value, y::AbstractExpr) = QolElemAtom(Constant(sqrt(x)), y)
  broadcast(::typeof(^), x::AbstractExpr, k::Int) = k==2 ? QolElemAtom(x, Constant(ones(size(x)))) : error("raising variables to powers other than 2 is not implemented")
end
