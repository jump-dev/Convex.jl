# Convex.jl - Convex Optimization in Julia

Convex.jl is a Julia package for [Disciplined Convex Programming](http://dcp.stanford.edu/) (DCP).

Convex.jl makes it easy to describe optimization problems in a natural,
mathematical syntax, and to solve those problems using a variety of different
(commercial and open-source) solvers.

Convex.jl can be used to solve:

 * linear programs
 * mixed-integer linear programs and mixed-integer second-order cone programs
 * DCP-compliant convex programs including
   * second-order cone programs (SOCP)
   * exponential cone programs
   * semidefinite programs (SDP)

## Resources for getting started

There are a few ways to get started with Convex:

 * Read the [Installation](@ref) guide
 * Read the introductory tutorial [Quick Tutorial](@ref)
 * Read the list of [Supported Operations](@ref)
 * Browse some of our examples

!!! tip
    Need help? Join the [community forum](https://jump.dev/forum)
    to search for answers to commonly asked questions.

    Before asking a question, make sure to read the post [make it easier to help you](https://discourse.julialang.org/t/psa-make-it-easier-to-help-you/14757),
    which contains a number of tips on how to ask a good question.

## How the documentation is structured

Having a high-level overview of how this documentation is structured will help
you know where to look for certain things.

* **Examples** contain worked examples of solving problems with Convex. Start
  here if you are new to Convex, or you have a particular problem class you want
  to model.

* The **Manual** contains short code-snippets that explain how to achieve
  specific tasks in Convex. Look here if you want to know how to achieve a
  particular task.

* The **Developer docs** section contains information for people contributing to
  Convex development. Don't worry about this section if you are using Convex to
  formulate and solve problems as a user.
