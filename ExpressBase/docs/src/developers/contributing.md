# [Contributing](@id contributing)

```@contents
Pages = ["contributing.md"]
```

Welcome! This document explains some ways you can contribute to `ExpressBase`.

## Code of conduct

This project and everyone participating in it is governed by the
["Contributor Covenant Code of Conduct"](https://github.com/MineralsCloud/.github/blob/main/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code.

## Join the community forum

First up, join the [community forum](https://github.com/MineralsCloud/ExpressBase.jl/discussions).

The forum is a good place to ask questions about how to use `ExpressBase`. You can also
use the forum to discuss possible feature requests and bugs before raising a
GitHub issue (more on this below).

Aside from asking questions, the easiest way you can contribute to `ExpressBase` is to
help answer questions on the forum!

## Improve the documentation

Chances are, if you asked (or answered) a question on the community forum, then
it is a sign that the [documentation](https://MineralsCloud.github.io/ExpressBase.jl/dev/) could be
improved. Moreover, since it is your question, you are probably the best-placed
person to improve it!

The docs are written in Markdown and are built using
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).
You can find the source of all the docs
[here](https://github.com/MineralsCloud/ExpressBase.jl/tree/master/docs).

If your change is small (like fixing typos, or one or two sentence corrections),
the easiest way to do this is via GitHub's online editor. (GitHub has
[help](https://help.github.com/articles/editing-files-in-another-user-s-repository/)
on how to do this.)

If your change is larger, or touches multiple files, you will need to make the
change locally and then use Git to submit a
[pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests).
(See [Contribute code to `ExpressBase`](@ref) below for more on this.)

## File a bug report

Another way to contribute to `ExpressBase` is to file
[bug reports](https://github.com/MineralsCloud/ExpressBase.jl/issues/new?template=bug_report.md).

Make sure you read the info in the box where you write the body of the issue
before posting. You can also find a copy of that info
[here](https://github.com/MineralsCloud/ExpressBase.jl/blob/master/.github/ISSUE_TEMPLATE/bug_report.md).

!!! tip
    If you're unsure whether you have a real bug, post on the
    [community forum](https://github.com/MineralsCloud/ExpressBase.jl/discussions)
    first. Someone will either help you fix the problem, or let you know the
    most appropriate place to open a bug report.

## Contribute code to `ExpressBase`

Finally, you can also contribute code to `ExpressBase`!

!!! warning
    If you do not have experience with Git, GitHub, and Julia development, the
    first steps can be a little daunting. However, there are lots of tutorials
    available online, including:
    * [GitHub](https://guides.github.com/activities/hello-world/)
    * [Git and GitHub](https://try.github.io/)
    * [Git](https://git-scm.com/book/en/v2)
    * [Julia package development](https://docs.julialang.org/en/v1/stdlib/Pkg/#Developing-packages-1)

Once you are familiar with Git and GitHub, the workflow for contributing code to
`ExpressBase` is similar to the following:

### Step 1: decide what to work on

The first step is to find an [open issue](https://github.com/MineralsCloud/ExpressBase.jl/issues)
(or open a new one) for the problem you want to solve. Then, _before_ spending
too much time on it, discuss what you are planning to do in the issue to see if
other contributors are fine with your proposed changes. Getting feedback early can
improve code quality, and avoid time spent writing code that does not get merged into
`ExpressBase`.

!!! tip
    At this point, remember to be patient and polite; you may get a _lot_ of
    comments on your issue! However, do not be afraid! Comments mean that people are
    willing to help you improve the code that you are contributing to `ExpressBase`.

### Step 2: fork `ExpressBase`

Go to [https://github.com/MineralsCloud/ExpressBase.jl](https://github.com/MineralsCloud/ExpressBase.jl)
and click the "Fork" button in the top-right corner. This will create a copy of
`ExpressBase` under your GitHub account.

### Step 3: install `ExpressBase` locally

Similar to [installation](@ref), open the Julia REPL and run:

```@repl
using Pkg
Pkg.update()
Pkg.develop("ExpressBase")
```

Then the package will be cloned to your local machine. On *nix systems, the default path is
`~/.julia/dev/ExpressBase` unless you modify the
[`JULIA_DEPOT_PATH`](http://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH-1)
environment variable. If you're on
Windows, this will be `C:\\Users\\<my_name>\\.julia\\dev\\ExpressBase`.
In the following text, we will call it `PKGROOT`.

Go to `PKGROOT`, start a new Julia session and run

```@repl
using Pkg
Pkg.instantiate()
```

to instantiate the project.

### Step 4: checkout a new branch

!!! note
    In the following, replace any instance of `GITHUB_ACCOUNT` with your GitHub
    username.

The next step is to checkout a development branch. In a terminal (or command
prompt on Windows), run:

```shell
cd ~/.julia/dev/ExpressBase

git remote add GITHUB_ACCOUNT https://github.com/GITHUB_ACCOUNT/ExpressBase.jl.git

git checkout master  # or main

git pull

git checkout -b my_new_branch
```

### Step 5: make changes

Now make any changes to the source code inside the `~/.julia/dev/ExpressBase`
directory.

Make sure you:

* Follow the [Style guide](@ref style) and [run `JuliaFormatter.jl`](@ref formatter)
* Add tests and documentation for any changes or new features

!!! tip
    When you change the source code, you'll need to restart Julia for the
    changes to take effect. This is a pain, so install
    [`Revise.jl`](https://github.com/timholy/Revise.jl).

### Step 6a: test your code changes

To test that your changes work, run the `ExpressBase` test-suite by opening Julia and
running:

```@repl
cd("~/.julia/dev/ExpressBase")
using Pkg
Pkg.activate(".")
Pkg.test()
```

!!! warning
    Running the tests might take a long time.

!!! tip
    If you're using `Revise.jl`, you can also run the tests by calling `include`:

    ```julia
    include("test/runtests.jl")
    ```

    This can be faster if you want to re-run the tests multiple times.

### Step 6b: test your documentation changes

Open Julia, then run:

```@repl
cd("~/.julia/dev/ExpressBase/docs")
using Pkg
Pkg.activate(".")
include("src/make.jl")
```

After a while, a folder `PKGROOT/docs/build` will appear. Open
`PKGROOT/docs/build/index.html` with your favorite browser, and have fun!

!!! warning
    Building the documentation might take a long time.

!!! tip
    If there's a problem with the tests that you don't know how to fix, don't
    worry. Continue to step 5, and one of the `ExpressBase` contributors will comment
    on your pull request telling you how to fix things.

### Step 7: make a pull request

Once you've made changes, you're ready to push the changes to GitHub. Run:

```shell
cd ~/.julia/dev/ExpressBase

git add .

git commit -m "A descriptive message of the changes"

git push -u GITHUB_ACCOUNT my_new_branch
```

Then go to [https://github.com/MineralsCloud/ExpressBase.jl/pulls](https://github.com/MineralsCloud/ExpressBase.jl/pulls)
and follow the instructions that pop up to open a pull request.

### Step 8: respond to comments

At this point, remember to be patient and polite; you may get a _lot_ of
comments on your pull request! However, do not be afraid! A lot of comments
means that people are willing to help you improve the code that you are
contributing to `ExpressBase`.

To respond to the comments, go back to step 5, make any changes, test the
changes in step 6, and then make a new commit in step 7. Your PR will
automatically update.

### Step 9: cleaning up

Once the PR is merged, clean-up your Git repository ready for the
next contribution!

```shell
cd ~/.julia/dev/ExpressBase

git checkout master

git pull
```

!!! note
    If you have suggestions to improve this guide, please make a pull request!
    It's particularly helpful if you do this after your first pull request
    because you'll know all the parts that could be explained better.

Thanks for contributing to `ExpressBase`!
