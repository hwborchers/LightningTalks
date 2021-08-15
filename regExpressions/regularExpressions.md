---
title: "Regular Expressions in R"
author: "Hans W. Borchers"
date: "September 2021"
output:
  html_document:
    keep_md: true
---



------------------------------------------------------------------------

### Introduction

Regular expression are patterns describing a certain amount of text.  
Their name comes from the mathematical theory on which they are based.

#### Resources

- [Regular Expressions Info](https://regular-expressions.info) Website  
  with *Tutorial*, Examples, Books, and much more

- [Regular Expressions](https://github.com/rstudio/cheatsheets/raw/master/regex.pdf) Cheatsheet contributed to RStudio

- Jeffrey Friedl. [Mastering Regular Expressions](https://www.google.de/books/edition/Mastering_Regular_Expressions/sshKXlr32-AC?hl=en&gbpv=1&dq=mastering+regular+expressions&printsec=frontcover).  
  3rd Edition, O'Reilly Media Inc., 2006

#### Caveats

* This talk is a 'compactification' of the regular expressions tutorial above.  
  There are many, many more introductions and tutorials on the Web.

* There are R functions treating regular expressions in Base R  
  and in the 'stringr' package. We will use the Base R version.

* The RE syntax is identical in almost all modern programming languages,  
  as all of them utilize the free and excellent [PCRE](https://www.pcre.org/)
  library providing  
  [PERL](https://perldoc.perl.org/perlre) compatible regular expressions.

#### An Example

    \s+[-+]?(?:[0-9]*[.])?[0-9]+(?:[eE][-+]?\d+)?\s+        What is this?

-----------------------------------------------------------------------

### Regular Expression syntax

#### Literal Characters

    # Literal characters:       a-z, A-Z, 0-9, _

    # Special characters:       \ ^ $ . | ? * + ( ) [ {

These special or 'metacharacters' have to be escaped with a backslash if used in a normal search.  
`] } ' "` are *not* special characters in normal text.

#### Character classes

A 'character class' is a set of characters that the regex machine is asked to match, for instance `[aeiou]` for all vowels.

Predefined/shorthand character classes are

    \w      words chars     [a-zA-Z0-9_]
    \d      digits          [0-9]
    \s      whitespace      [ \t\r\f\n]

Negated character classes are introduced with `^`, that is `[^aeiou]` matches all non-vocals. Predefined negated shorthand classes are

    \W      [^\w]
    \D      [^\d]
    \S      [^\s]

Character classes can be combined into, e.g., `[\s\d]`. But note that `[\D\S]` is not the same as `[^\d\s]`.

#### Dots and anchors

    .       matches *any* character (except newline `\n`)

    ^       begin, position before the first character
    $       end, right after the last character

    \b      word boundary (incl. string start and end)
    \B      negated \b

#### Alternation

    |       "(cat|dog)" matches cat or dog

#### Optional characters

    ?       Nov(ember)? matches Nov and November -- greedy !
    ??      Nov(ember)?? matches Nov even if November is present  - lazy !

#### Repetition

    ?       zero or once        = {0,1}
    *       zero or more        = {0,}
    +       once or more        = {1,}
    {m}     exactly m times
    {m,}    m or more times
    {m,n}   betweeen m and n times ({,n} not valid)

Watch out for *greediness*; the regex engine tries to find a maximal match.
Repetition can be forced to be *lazy* by appending a '?'. To find a HTML tag, `re = r"<.+>"` is wrong, use `r"<.+?>"` or  better `r"<[^>]+>"`.

#### Grouping and capturing

    (...)               captured groups
    (?:...)             non-capturing groups

    (?<name>...)        named capture groups

#### Backreferences

Backreferences match the same text as previously matched by a capturing group.

    \1, \2, \3, ...     matches the first, second, ... capturing group

    \k<name>            (or \k'name') matches named backreferences

#### Lookahead and lookbehind

    q(?!u)      negative lookahead: q *not* followed by u
    q(?=u)      positive lookahead: q followed by u,
                                    but u is not part of the match

    (?<!a)b     negative lookbehind:  matches b not preceded by an a
    (?<=a)b     positive lookbehind:  b preceded by an a,
                                      but a is not part of the match

    a\Kb        simplified lookbehind:  keeps text matched so far
                (only Perl, PCRE)       out of the overall regex match

#### Specifying modes

Start a regular expression with one or many of these options:

    (?i)        makes regex case-insensitive
    (?x)        turns on free-spacing mode
    (?xx)       free-spacing mode in character classes
    (?s)        single-line mode ('.' matches '\n')
    (?m)        multi-line mode ('^' and '$' match on each line)
    (?n)        turns unnamed grops non-capturing
    (?U)        turns on ungreedy mode            (PCRE only)
    (?X)        makes '\x' an error if not valid  (PCRE only)

------------------------------------------------------------------------

### Functions for regular expressions

#### Base R

R does have a "regular expression" string class (like, e.g., Python or Julia), a backslash '\' must be doubled (Ex. "\\w") to survive, but not in r"(\w)". Also, r"{...}" or r"[...]" can be used.

**grep and grepl**

    grep[l](pat, str, ignore.case = FALSE,
            perl = FALSE, value = FALSE, fixed = FALSE)

**regexpr and gregexpr**

    [g]regexpr(pat, str, ignore.case = FALSE, perl = FALSE,
               fixed = FALSE, useBytes = FALSE)

**strsplit**

    strsplit(str, split, fixed = FALSE, perl = FALSE, useBytes = FALSE)

**sub and gsub**

    [g]sub(pat, replace, str, ignore.case = FALSE,
           perl = FALSE, fixed = FALSE, useBytes = FALSE)

**regmatches**

    m <- regexpr(pat, str)
    regmatches(str, m, invert = FALSE)
    regmatches(str, m, invert = FALSE) <- value

------------------------------------------------------------------------

### useful Examples

    ^\s*  \s*^


------------------------------------------------------------------------

### More R packages

#### The 'stringr' package

    str_detect(str, pat)
    str_locate()                str_locate_all()
    str_extract()               str_extract_all()
    str_match()                 str_match_all()
    str_replace(str, pat, rep)  str_replace_all
    str_split()

#### re2

're2' provides an interface to the Google regular-expression library.

"PCRE relies on recursive backtracking. For this approach time complexity can be exponential. In contrast, re2 uses finite automata-based techniques for regular expression matching, guaranteeing linear time execution and a fixed stack footprint."

#### regextestR

`run_app()` starts a 'Shiny' application "R Regex Tester" with input fields for a matching pattern and a test string. Successful matches will be highlighted.

#### namedCapture

"User-friendly wrappers for named capture regular expressions." With extensive vignettes, and covered by an R Journal article in 2019.

#### RVerbalExpressions

"Create regular expressions easily ... using grammar and functionality inspired by [VerbalExpressions](https://github.com/VerbalExpressions). Usage of the 'pipe' is encouraged to build expressions in a chain-like fashion."

#### rematch and rematch2

'rematch2' wraps "on 'regexpr' and 'gregexpr' to return the match results in tidy data frames."

#### rex and rebus

'rex' provides "a user-friendly interface for the construction of regular expressions."

'rebus' "builds regular expressions piece by piece using human readable code."
