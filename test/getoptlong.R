library(GetoptLong)

# cutoff = 0.05
# verbose = TRUE
# GetoptLong(
#   "number=i", "Number of items.",
#   "cutoff=f", "Cutoff for filtering results.",
#   "verbose!", "Print message."
# )
# 
# 
# print(number)

count = 1
number = 0.1
array = c(1, 2)
hash = list("foo" = "a", "bar" = "b")
verbose = TRUE
GetoptLong.options(help_style = "two-column")
GetoptLong(
  ...
)
