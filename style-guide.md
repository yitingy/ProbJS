# Style guide

## Problematic constructs

Due to lack of stack extension:
- for
- ?list-comprehension
- ...


## Mutation

1. While mutation is not a problem for lightweight MH (as in Bher), the long-range dependencies introduced by mutation can cause for more sophisticated MH algorithms.

2. In thinking about the language, consider the question: in the proposed first book (which is just forward sampling), how are we going beyond what the procedural community has been doing all along?  "Programs with random choices" seems to be what they have been doing before. One thing we could be adding is a clean interpretation--if subexprs (+ environments) have their own marginal distributions, they factor from the rest of the program and one can think about replacing them with other distributions. However, this is only the case if programs are functional.
