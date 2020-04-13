# OPERA: Ordering Poset Elements by Recursive Amalgamation

Risk stratifying patients using established risk factors is a perpetual research theme. 

Statistically, this task is a supervised clustering or stratification problem where the outcome is usually time to event (such as recurrence or death) that is subject to possible censoring due to loss to follow up. But the ordinal structure of staging variables often arise the paritial ordinal problem in the staging problem. 
For example, the interaction terms of two ordinal categorical variables are patial ordered, with comparable and incomparable pairs. 

we generalize the partially order relationship appeared in such A $\times$ B interaction to all partially ordered set (**poset**), which enable us to solve the similar staging problems with prior knowledge on ordering. We propose a recursive algorithm to order and amalgamate the poset elements into relatively few strata. This method can generate more flexible splitting patterns of the categories while preserving the partial order relationship of the poset.

## Installation

library(devtools)

install_github("WangTJ/opera")

## Structure of Repositories

For your inference, the functions help documents can be found in <https://github.com/WangTJ/opera/tree/master/man>, the example code can be found in <https://github.com/WangTJ/opera/tree/master/analysis/example.R>
