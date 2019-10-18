'-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?'
One piece at a time:

-?        optionally matches a negative sign (zero or one negative signs)
\ *       matches any number of spaces (to allow for formatting variations like - 2.3 or -2.3)
[0-9]+    matches one or more digits
\.?       optionally matches a period (zero or one periods)
[0-9]*    matches any number of digits, including zero
(?: ... ) groups an expression, but without forming a "capturing group" (look it up)
[Ee]      matches either "e" or "E"
\ *       matches any number of spaces (to allow for formats like 2.3E5 or 2.3E 5)
-?        optionally matches a negative sign
\ *       matches any number of spaces
[0-9]+    matches one or more digits
?         makes the entire non-capturing group optional (to allow for the presence or absence of the exponent - 3000 or 3E3
