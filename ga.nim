assert NimMajor >= 2, "must be nim 2.0"

# imports

import malebolgia
import std/[math, random, terminal, strformat, algorithm, sequtils, sugar]

# init

randomize()

# functions related for sampling mutation length from power law distribution

proc cumsum(x: seq[float]): seq[float] =
  var sum = 0.0

  for el in x:
    sum += el
    result.add(sum)

proc sample(cdf: seq[float]): int =
  var n = cdf.len
  var x = rand(1.0)

  if x < cdf[0]:
    return 1

  var lo = 0
  var hi = n - 1

  while (hi - lo) > 1:
    var mid = floor_div(lo + hi, 2)

    if x > cdf[mid]:
      lo = mid
    else:
      hi = mid

  return hi + 1

proc log1p(x: float): float =
  return ln(1 + x)

proc sample_waiting_time(prob: float): int =
  var rand_prob = rand(1.0)
  return 1 + (log1p(-rand_prob) / log1p(-prob)).floor.int

# gene

type
  Gene = object
    code: string
    cost: int = 1e7.int
  GeneRef = ref Gene

proc init_random_code(gene: GeneRef, length: int) =
  assert length > 2

  while len(gene.code) < length:
    gene.code = gene.code & rand(255).chr

# finds a random index in the code and adds or subtracts 1 from the ord randomly

type
  MutateOperatorType = enum
    DefaultMutateOperator,
    FastMutateOperator

  MutateOperator = object
    gene_length: int = 0

    case kind: MutateOperatorType
      of DefaultMutateOperator:
        mutate_rate: float = 0.1

      of FastMutateOperator:
        beta: float = 1.5
        power_law: seq[float]

  MutateOperatorRef = ref MutateOperator

proc init(mutate_operator: MutateOperatorRef, gene_length: int) =
  mutate_operator.gene_length = gene_length

  case mutate_operator.kind:
  of FastMutateOperator:
    var beta = mutate_operator.beta
    var half_gene_length = floor_div(gene_length, 2)

    var power_law = to_seq(1..half_gene_length).map(x => x.float)

    power_law = power_law.map(value => value.pow(-beta))
    power_law = cumsum(power_law)
    power_law = power_law.map(value => value / power_law[^1])

    mutate_operator.power_law = power_law

  else:
    discard

proc mutate(mutate_operator: MutateOperatorRef, gene: GeneRef) =

  var code_seq_char = cast[seq[char]](gene.code)
  var gene_length = gene.code.len

  case mutate_operator.kind:
  of DefaultMutateOperator:

    for i in 0..<gene_length:
      if rand(1.0) > mutate_operator.mutate_rate:
        continue

      let mutation = (if rand(1.0) < 0.5: 1 else: -1)

      let char_at_index = code_seq_char[i]

      var char_code = char_at_index.ord + mutation
      char_code = char_code.max(0).min(255)

      code_seq_char[i] = char_code.chr

  else:

    let mutate_rate = sample(mutate_operator.power_law) / gene_length
    var index = sample_waiting_time(mutate_rate) - 1

    while index < gene_length:

      let mutation = (if rand(1.0) < 0.5: 1 else: -1)
      let char_at_index = code_seq_char[index]
      var char_code = char_at_index.ord + mutation

      char_code = char_code.max(0).min(255)
      code_seq_char[index] = char_code.chr

      index += sample_waiting_time(mutate_rate)

  gene.code = cast[string](code_seq_char)


# recombine genetic codes

proc mate(gene1, gene2: GeneRef): GeneRef = 
  assert gene1.code.len == gene2.code.len
  let pivot = rand(1..(gene1.code.len - 2))
  let new_code = gene1.code[0..<pivot] & gene2.code[pivot..^1]
  return GeneRef(code: new_code)

func calc_cost(code: string, goal: string): int =
  assert len(code) == len(goal)

  var cost: int = 0

  for i in 0..<len(goal):
    let diff: int = goal[i].ord - code[i].ord
    cost += pow(diff.float, 2).int

  return cost

proc cmp_genes(gene1, gene2: GeneRef): int =
  gene1.cost - gene2.cost

# population

type
  Population = object
    size: int
    pool: seq[GeneRef]
    goal: string 
    generation: int
    keep_fittest_frac: float = 1.0
    solved: bool = false
    mutate_prob: float = 0.5
    mutate_operator: MutateOperatorRef

  PopulationRef = ref Population

proc sort_by_fitness(pop: PopulationRef) =
  pop.pool.sort(cmp_genes, Ascending)

proc init_population(pop: PopulationRef) =

  var gene_length = pop.goal.len

  for _ in 0..<pop.size:
    let gene = GeneRef()
    gene.init_random_code(gene_length)
    pop.pool.add(gene)

  sort_by_fitness(pop)
  pop.mutate_operator.init(gene_length)

proc keep_fittest(pop: PopulationRef) =
  pop.sort_by_fitness()
  let num_survive: int = (pop.size.float * pop.keep_fittest_frac).float.int
  pop.pool = pop.pool[0..<num_survive]

proc breed_fittest(pop: PopulationRef) =
  assert pop.pool.len >= 2

  while len(pop.pool) < pop.size:
    pop.pool.shuffle()
    let child = mate(pop.pool[0], pop.pool[1])
    assert child.code.len == pop.goal.len
    pop.pool.add(child)

# next generation

proc next_generation(pop: PopulationRef) =

  if pop.solved:
    return

  # if population not started, initialize with random genes

  if pop.generation == 0:
    init_population(pop)

  # keep only the fittest

  pop.keep_fittest()

  # breed the fittest

  pop.breed_fittest()

  # mutate

  for gene in pop.pool:
    if rand(1.0) < pop.mutate_prob:
      pop.mutate_operator.mutate(gene)

  # calculate fitness

  var master = create_master()
  var costs = new_seq[int](pop.pool.len)

  master.await_all:
    for i in 0..<pop.pool.len:
      master.spawn calc_cost(pop.pool[i].code, pop.goal) -> costs[i]

  for i in 0..<pop.pool.len:
    let gene = pop.pool[i]
    gene.cost = costs[i]

  # sort

  pop.sort_by_fitness()

  # determine solved

  pop.solved = pop.pool[0].cost == 0

  # increment generation

  pop.generation += 1

# display population of genes

proc display(pop: PopulationRef) =
  erase_screen()

  echo &"generation: {pop.generation}", "\n"

  for gene in pop.pool:
    echo &"{gene.code} ({gene.cost})"

  echo "\n"

proc main() = 

  let trials = 100
  let use_fast_mutate = true
  var avg_generation = 0.0

  var mutate_operator: MutateOperatorRef

  if use_fast_mutate:
    mutate_operator = MutateOperatorRef(
      kind: FastMutateOperator,
      beta: 1.2
    )
  else:
    mutate_operator = MutateOperatorRef(
      kind: DefaultMutateOperator,
      mutate_rate: 0.1
    )

  # do many trials and average generation at solving

  for trial in 1..trials:

    # create a gene population

    let population = PopulationRef(
      goal: "Attention is all you need",
      size: 25,
      keep_fittest_frac: 0.25,
      mutate_prob: 0.9,
      mutate_operator: mutate_operator
    )

    # while not solved, do another generation

    while not population.solved:
      population.next_generation()

    avg_generation += population.generation / trials

    echo &"trial {trial} completed at generation {population.generation}"

  echo &"average generation: {avg_generation.int}"

when is_main_module:
  main()
