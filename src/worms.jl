const WormIndex = Int
const StateIndex= UInt8

worm_action(worm::Integer, state::Integer, basis_size::Integer) =
    (state + worm - 1) % basis_size + 1
worm_inverse(worm::Integer, basis_size::Integer) = basis_size - worm
worm_count(basis_size::Integer) = basis_size - 1

# worm_action(worm::Integer, state::Integer, ::Integer) = xor(state - 1, worm) + 1
# worm_inverse(worm::Integer, ::Integer) = worm
# worm_count(basis_size::Integer) = basis_size - 1
