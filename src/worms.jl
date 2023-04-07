const WormIdx = Int
const StateIdx = UInt8

worm_action(worm::WormIdx, state::StateIdx, ::Integer) = ((state - 1)^(worm)) + 1
worm_inverse(worm::WormIdx, ::Integer) = worm
worm_count(basis_size::Integer) = basis_size - 1
