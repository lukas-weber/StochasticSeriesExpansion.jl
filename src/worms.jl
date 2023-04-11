const WormIdx = Int
const StateIdx = UInt8

worm_action(worm::Integer, state::Integer, ::Integer) = xor(state - 1, worm) + 1
worm_inverse(worm::Integer, ::Integer) = worm
worm_count(basis_size::Integer) = basis_size - 1
