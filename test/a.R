x <- rnorm(100, 0, 1)
y <- runif(100, 0, 1)
?rnorm


fun <- function(var1, var2) {
  var0 <- var1 + var2
  var0
}

lm(y ~ x)
?plot
?plot
# section 1

# section --------------------------------------------------------------------
library(ggplot2)
ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_point()

mtcars

d <- mtcars

View(mtcars)

library(plotly)
ggpenguins <- qplot(
  bill_length_mm,
  body_mass_g,
  data = palmerpenguins::penguins,
  color = species,
)
ggplotly(ggpenguins)
DT::datatable(mtcars)
a <- TRUE
b <- FALSE
library(lintr)
lint(
  text = "mean(x, trim = 0.2, na.rm = TRUE)",
  linters = commas_linter()
)

paletteer::paletteer_c("viridis::viridis", n = 10)
