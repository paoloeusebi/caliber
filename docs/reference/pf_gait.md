# NHANES 2001-2002 gait speed and mobility items, ages 60+

Gait speed from the timed 20-ft walk and seven lower-body mobility items
of the physical functioning questionnaire for NHANES 2001-2002
participants aged 60+ with a valid walk time and at least one item
answered. Analyses are unweighted. \[gait_data()\] returns these data,
optionally filtered or recoded, together with the item matrix.

## Usage

``` r
pf_gait
```

## Format

A data frame with 1426 rows and 15 columns:

- SEQN:

  NHANES respondent sequence number

- age:

  age in years at screening (60 to 85; 85 = 85 or older)

- sex:

  1 = male, 2 = female

- walk_time:

  time to walk 20 ft, seconds

- gait_speed:

  gait speed, m/s (6.096 m / walk_time)

- walking_aid:

  1 = walking-aid user (needs special equipment to walk, used a device
  in the walk test, or could not walk unassisted), 0 = not

- C1:

  difficulty walking a quarter mile (0 = no difficulty, 1 = some, 2 =
  much, 3 = unable; NA if missing or not asked)

- C2:

  difficulty walking up 10 steps (coded as C1)

- C3:

  difficulty walking between rooms on the same floor (coded as C1)

- C4:

  difficulty standing up from an armless chair (coded as C1)

- C5:

  difficulty getting in and out of bed (coded as C1)

- C6:

  difficulty stooping, crouching or kneeling (coded as C1)

- C7:

  difficulty standing for long periods (coded as C1)

- n_items_answered:

  number of items C1 to C7 answered (1 to 7)

- needs_equipment:

  1 = needs special equipment to walk (PFQ055 = 1); C1 and C2 are not
  asked of these people

## Source

NHANES 2001-2002 public-use files MSX_B (lower extremity functioning),
PFQ_B (physical functioning) and DEMO_B (demographics), National Center
for Health Statistics, \<https://wwwn.cdc.gov/nchs/nhanes/\>. Built by
\`inst/scripts/gait_data.R\`.
