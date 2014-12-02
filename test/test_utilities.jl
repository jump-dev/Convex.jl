using Base.Test
using Convex

# length and size
x = Variable(2,3)
@test length(x)==6
@test size(x)==(2,3)
@test size(x,1)==2
@test size(x,2)==3

x = Variable(3)
@test length(x)==3
@test size(x)==(3,1)

x = Variable()
@test length(x)==1
@test size(x)==(1,1)