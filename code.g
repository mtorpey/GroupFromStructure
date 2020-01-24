NonTrivialSemidirectProduct := function(N, H)
  local homos, kernel_sizes, pos, semidir;
  homos := AllHomomorphismClasses(H, AutomorphismGroup(N));;  # all homos
  kernel_sizes := List(homos, homo -> Size(Kernel(homo)));
  pos := PositionProperty(kernel_sizes, s -> s < Size(H));
  if pos = fail then
    return fail;  # only trivial/universal homomorphisms exist
  fi;
  semidir := SemidirectProduct(H, homos[pos], N);  # create product
  return semidir;
end;

GroupByStructureDescription := function(str)
  local from, depth, ind, G1, G2, N, H, semidir, args, d, q, n;
  # Remove spaces
  str := Filtered(str, c -> c <> ' ');
  
  # Handle any parentheses at the start
  from := 0;
  if str[1] = '(' then
    depth := 1;
    from := 1;
    repeat
      if from >= Length(str) then
        ErrorNoReturn("Mismatched parentheses: ", str, "\n");
      fi;
      if str[from + 1] = '(' then
        depth := depth + 1;
      elif str[from + 1] = ')' then
        depth := depth - 1;
      fi;
      from := from + 1;
    until depth = 0;
    if from >= Length(str) then
      # parentheses around whole statement
      str := str{[2 .. Length(str) - 1]};
      from := 0;
    fi;
  fi;
  
  # Handle direct products
  ind := Position(str, 'x', from);
  if ind <> fail then
    G1 := GroupByStructureDescription(str{[1 .. ind - 1]});
    G2 := GroupByStructureDescription(str{[ind + 1 .. Length(str)]});
    return DirectProduct(G1, G2);
  fi;
  
  # Handle semidirect products
  ind := Position(str, ':', from);
  if ind <> fail then
    N := GroupByStructureDescription(str{[1 .. ind - 1]});
    H := GroupByStructureDescription(str{[ind + 1 .. Length(str)]});
    semidir := NonTrivialSemidirectProduct(N, H);  # try to get a non-direct one
    if semidir = fail then
      return DirectProduct(N, H);  # direct products are semidirect products!
    fi;
    return semidir;
  fi;
  
  # Handle dots

  # Handle powers (short only)

  # Atom
  if str = "QD16" then
    return Group([(1,2,4,7)(3,5,8,6), (2,5)(3,8)(6,7),
                  (1,3,4,8)(2,6,7,5), (1,4)(2,7)(3,8)(5,6)]);
  elif Length(str) >= 3 and str{[1 .. 3]} = "SL(" and str[Length(str)] = ')' then
    args := SplitString(str{[4 .. Length(str) - 1]}, ",");
    d := Int(args[1]);
    q := Int(args[2]);
    return SpecialLinearGroup(d, q);
  elif str[1] = 'A' then
    n := Int(str{[2..Length(str)]});
    return AlternatingGroup(n);
  elif str[1] = 'C' then
    n := Int(str{[2..Length(str)]});
    return CyclicGroup(IsPermGroup, n);
  elif str[1] = 'D' then
    n := Int(str{[2..Length(str)]});
    return DihedralGroup(IsPermGroup, n);
  elif str[1] = 'Q' then
    n := Int(str{[2..Length(str)]});
    return QuaternionGroup(IsPermGroup, n);
  elif str[1] = 'S' then
    n := Int(str{[2..Length(str)]});
    return SymmetricGroup(n);
  elif str[1] in "0123456789" then
    n := Int(str);
    return CyclicGroup(IsPermGroup, n);
  fi;
  Error("Couldn't match ", str);
end;
