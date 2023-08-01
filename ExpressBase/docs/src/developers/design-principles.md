# Design Principles

```@contents
Pages = ["design-principles.md"]
Depth = 2
```

We adopt some [`SciML`](https://sciml.ai/) design [guidelines](https://github.com/SciML/SciMLStyle)
here. Please read it before contributing!

## Consistency vs adherence

According to PEP8:

> A style guide is about consistency. Consistency with this style guide is important.
> Consistency within a project is more important. Consistency within one module or function is the most important.
>
> However, know when to be inconsistent—sometimes style guide recommendations just aren't
> applicable. When in doubt, use your best judgment. Look at other examples and decide what
> looks best. And don’t hesitate to ask!

## Community contribution guidelines

For a comprehensive set of community contribution guidelines, refer to [ColPrac](https://github.com/SciML/ColPrac).
A relevant point to highlight PRs should do one thing. In the context of style, this means that PRs which update
the style of a package's code should not be mixed with fundamental code contributions. This separation makes it
easier to ensure that large style improvement are isolated from substantive (and potentially breaking) code changes.

## Open source contributions are allowed to start small and grow over time

If the standard for code contributions is that every PR needs to support every possible input type that anyone can
think of, the barrier would be too high for newcomers. Instead, the principle is to be as correct as possible to
begin with, and grow the generic support over time. All recommended functionality should be tested, any known
generality issues should be documented in an issue (and with a `@test_broken` test when possible).

## Generic code is preferred unless code is known to be specific

For example, the code:

```@repl
function f(A, B)
    for i in 1:length(A)
        A[i] = A[i] + B[i]
    end
end
```

would not be preferred for two reasons. One is that it assumes `A` uses one-based indexing, which would fail in cases
like [OffsetArrays](https://github.com/JuliaArrays/OffsetArrays.jl) and [FFTViews](https://github.com/JuliaArrays/FFTViews.jl).
Another issue is that it requires indexing, while not all array types support indexing (for example,
[CuArrays](https://github.com/JuliaGPU/CuArrays.jl)). A more generic compatible implementation of this function would be
to use broadcast, for example:

```@repl
function f(A, B)
    @. A = A + B
end
```

which would allow support for a wider variety of array types.

## Internal types should match the types used by users when possible

If `f(A)` takes the input of some collections and computes an output from those collections, then it should be
expected that if the user gives `A` as an `Array`, the computation should be done via `Array`s. If `A` was a
`CuArray`, then it should be expected that the computation should be internally done using a `CuArray` (or appropriately
error if not supported). For these reasons, constructing arrays via generic methods, like `similar(A)`, is preferred when
writing `f` instead of using non-generic constructors like `Array(undef,size(A))` unless the function is documented as
being non-generic.

## Trait definition and adherence to generic interface is preferred when possible

Julia provides many interfaces, for example:

- [Iteration](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration)
- [Indexing](https://docs.julialang.org/en/v1/manual/interfaces/#Indexing)
- [Broadcast](https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting)

Those interfaces should be followed when possible. For example, when defining broadcast overloads,
one should implement a `BroadcastStyle` as suggested by the documentation instead of simply attempting
to bypass the broadcast system via `copyto!` overloads.

When interface functions are missing, these should be added to Base Julia or an interface package,
like [ArrayInterface.jl](https://github.com/JuliaArrays/ArrayInterface.jl). Such traits should be
declared and used when appropriate. For example, if a line of code requires mutation, the trait
`ArrayInterface.ismutable(A)` should be checked before attempting to mutate, and informative error
messages should be written to capture the immutable case (or, an alternative code which does not
mutate should be given).

One example of this principle is demonstrated in the generation of Jacobian matrices. In many scientific
applications, one may wish to generate a Jacobian cache from the user's input `u0`. A naive way to generate
this Jacobian is `J = similar(u0,length(u0),length(u0))`. However, this will generate a Jacobian `J` such
that `J isa Matrix`.

## Macros should be limited and only be used for syntactic sugar

Macros define new syntax, and for this reason they tend to be less composable than other coding styles
and require prior familiarity to be easily understood. One principle to keep in mind is, "can the person
reading the code easily picture what code is being generated?". For example, a user of Soss.jl may not know
what code is being generated by:

```julia
@model (x, α) begin
    σ ~ Exponential()
    β ~ Normal()
    y ~ For(x) do xj
        Normal(α + β * xj, σ)
    end
    return y
end
```

and thus using such a macro as the interface is not preferred when possible. However, a macro like
[`@muladd`](https://github.com/SciML/MuladdMacro.jl) is trivial to picture on a code (it recursively
transforms `a*b + c` to `muladd(a,b,c)` for more
[accuracy and efficiency](https://en.wikipedia.org/wiki/Multiply%E2%80%93accumulate_operation)), so using
such a macro for example:

```julia
julia> @macroexpand(@muladd k3 = f(t + c3 * dt, @. uprev + dt * (a031 * k1 + a032 * k2)))
:(k3 = f((muladd)(c3, dt, t), (muladd).(dt, (muladd).(a032, k2, (*).(a031, k1)), uprev)))
```

is recommended. Some macros in this category are:

- `@inbounds`
- [`@muladd`](https://github.com/SciML/MuladdMacro.jl)
- `@view`
- [`@named`](https://github.com/SciML/ModelingToolkit.jl)
- `@.`
- [`@..`](https://github.com/YingboMa/FastBroadcast.jl)

Some performance macros, like `@simd`, `@threads`, or
[`@turbo` from LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl),
make an exception in that their generated code may be foreign to many users. However, they still are
classified as appropriate uses as they are syntactic sugar since they do (or should) not change the behavior
of the program in measurable ways other than performance.

## Errors should be caught as high as possible, and error messages should be contextualized for newcomers

Whenever possible, defensive programming should be used to check for potential errors before they are encountered
deeper within a package. For example, if one knows that `f(u0,p)` will error unless `u0` is the size of `p`, this
should be caught at the start of the function to throw a domain specific error, for example "parameters and initial
condition should be the same size".

## Subpackaging and interface packages is preferred over conditional modules via Requires.jl

Requires.jl should be avoided at all costs. If an interface package exists, such as
[ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl) for defining automatic differentiation
rules without requiring a dependency on the whole ChainRules.jl system, or
[RecipesBase.jl](https://github.com/JuliaPlots/RecipesBase.jl) which allows for defining Plots.jl
plot recipes without a dependency on Plots.jl, a direct dependency on these interface packages is
preferred.

Otherwise, instead of resorting to a conditional dependency using Requires.jl, it is
preferred one creates subpackages, i.e. smaller independent packages kept within the same Github repository
with independent versioning and package management. An example of this is seen in
[Optimization.jl](https://github.com/SciML/Optimization.jl) which has subpackages like
[OptimizationBBO.jl](https://github.com/SciML/Optimization.jl/tree/master/lib/OptimizationBBO) for
BlackBoxOptim.jl support.

Some important interface packages to know about are:

- [ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl)
- [RecipesBase.jl](https://github.com/JuliaPlots/RecipesBase.jl)
- [ArrayInterface.jl](https://github.com/JuliaArrays/ArrayInterface.jl)
- [CommonSolve.jl](https://github.com/SciML/CommonSolve.jl)
- [SciMLBase.jl](https://github.com/SciML/SciMLBase.jl)

## Functions should either attempt to be non-allocating and reuse caches, or treat inputs as immutable

Mutating codes and non-mutating codes fall into different worlds. When a code is fully immutable,
the compiler can better reason about dependencies, optimize the code, and check for correctness.
However, many times a code making the fullest use of mutation can outperform even what the best compilers
of today can generate. That said, the worst of all worlds is when code mixes mutation with non-mutating
code. Not only is this a mishmash of coding styles, it has the potential non-locality and compiler
proof issues of mutating code while not fully benefiting from the mutation.

## Out-of-place and immutability is preferred when sufficient performant

Mutation is used to get more performance by decreasing the amount of heap allocations. However,
if it's not helpful for heap allocations in a given spot, do not use mutation. Mutation is scary
and should be avoided unless it gives an immediate benefit. For example, if
matrices are sufficiently large, then `A*B` is as fast as `mul!(C,A,B)`, and thus writing
`A*B` is preferred (unless the rest of the function is being careful about being fully non-allocating,
in which case this should be `mul!` for consistency).

Similarly, when defining types, using `struct` is preferred to `mutable struct` unless mutating
the struct is a common occurrence. Even if mutating the struct is a common occurrence, see whether
using [Setfield.jl](https://github.com/jw3126/Setfield.jl) is sufficient. The compiler will optimize
the construction of immutable structs, and thus this can be more efficient if it's not too much of a
code hassle.

## Tests should attempt to cover a wide gamut of input types

Code coverage numbers are meaningless if one does not consider the input types. For example, one can
hit all the code with `Array`, but that does not test whether `CuArray` is compatible! Thus, it's
always good to think of coverage not in terms of lines of code but in terms of type coverage. A good
list of number types to think about are:

- `Float64`
- `Float32`
- `Complex`
- [`Dual`](https://github.com/JuliaDiff/ForwardDiff.jl)
- `BigFloat`

Array types to think about testing are:

- `Array`
- [`OffsetArray`](https://github.com/JuliaArrays/OffsetArrays.jl)
- [`CuArray`](https://github.com/JuliaGPU/CUDA.jl)

## When in doubt, a submodule should become a subpackage or separate package

Keep packages to one core idea. If there's something separate enough to be a submodule, could it
instead be a separate well-tested and documented package to be used by other packages? Most likely
yes.

## Globals should be avoided whenever possible

Global variables should be avoided whenever possible. When required, global variables should be
constants and have an all uppercase name separated with underscores (e.g. `MY_CONSTANT`). They should be
defined at the top of the file, immediately after imports and exports but before an `__init__` function.
If you truly want mutable global style behavior you may want to look into mutable containers.

## Type-stable and type-grounded code is preferred wherever possible

Type-stable and type-grounded code helps the compiler create not only more optimized code, but also
faster to compile code. Always keep containers well-typed, functions specializing on the appropriate
arguments, and types concrete.

## Closures should be avoided whenever possible

Closures can cause accidental type instabilities that are difficult to track down and debug; in the
long run it saves time to always program defensively and avoid writing closures in the first place,
even when a particular closure would not have been problematic. A similar argument applies to reading
code with closures; if someone is looking for type instabilities, this is faster to do when code does
not contain closures.
See examples [here](https://discourse.julialang.org/t/are-closures-should-be-avoided-whenever-possible-still-valid-in-julia-v1-9/95893/5).

Furthermore, if you want to update variables in an outer scope, do so explicitly with `Ref`s or self
defined structs.

## Numerical functionality should use the appropriate generic numerical interfaces

While you can use `A\b` to do a linear solve inside a package, that does not mean that you should.
This interface is only sufficient for performing factorizations, and so that limits the scaling
choices, the types of `A` that can be supported, etc. Instead, linear solves within packages should
use LinearSolve.jl. Similarly, nonlinear solves should use NonlinearSolve.jl. Optimization should use
Optimization.jl. Etc. This allows the full generic choice to be given to the user without depending
on every solver package (effectively recreating the generic interfaces within each package).

## Functions should capture one underlying principle

Functions mean one thing. Every dispatch of `+` should be "the meaning of addition on these types".
While in theory you could add dispatches to `+` that mean something different, that will fail in
generic code for which `+` means addition. Thus, for generic code to work, code needs to adhere to
one meaning for each function. Every dispatch should be an instantiation of that meaning.

## Internal choices should be exposed as options whenever possible

Whenever possible, numerical values and choices within scripts should be exposed as options
to the user. This promotes code reusability beyond the few cases the author may have expected.

## Prefer code reuse over rewrites whenever possible

If a package has a function you need, use the package. Add a dependency if you need to. If the
function is missing a feature, prefer to add that feature to said package and then add it as a
dependency. If the dependency is potentially troublesome, for example because it has a high
load time, prefer to spend time helping said package fix these issues and add the dependency.
Only when it does not seem possible to make the package "good enough" should using the package
be abandoned. If it is abandoned, consider building a new package for this functionality as you
need it, and then make it a dependency.

## Prefer to not shadow functions

Two functions can have the same name in Julia by having different namespaces. For example,
`X.f` and `Y.f` can be two different functions, with different dispatches, but the same name.
This should be avoided whenever possible. Instead of creating `MyPackage.sort`, consider
adding dispatches to `Base.sort` for your types if these new dispatches match the underlying
principle of the function. If it doesn't, prefer to use a different name. While using `MyPackage.sort`
is not conflicting, it is going to be confusing for most people unfamiliar with your code,
so `MyPackage.special_sort` would be more helpful to newcomers reading the code.
