from typing import List, Optional
import math

def ageTransfer(
  age: float,
  adultAge: Optional[int] = 20
):
  age = (age+1)/(1+adultAge)
  if age <= 1:
    age = math.log(age)
  else:
    age = age - 1
  return age


def antiAgeTransfer(
     age: float,
     adultAge: Optional[int] = 20
  ):
    if age < 0:
        age = (1 + adultAge)*math.exp(age) - 1
    else:
        age = (1 + adultAge)*age + adultAge
    return age