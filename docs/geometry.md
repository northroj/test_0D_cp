

# [geometry]

Defines a cartesian mesh for the geometry (just one large x,y,z grid).

| Field | Type | Description |
|--------|------|-------------|
| `x` | list[float] or float float int mode | Defines x binning (increasing bin bounds OR start end N lin/log) |
| `y` | list[float] or float float int mode | Defines y binning (increasing bin bounds OR start end N lin/log) |
| `z` | list[float] or float float int mode | Defines z binning (increasing bin bounds OR start end N lin/log) |
| `fill` | int | Fill every cell with mat_id ? |
| `block` | 7 ints | Fills a block with a given material [-x_idx, +x_idx, -y_idx, +y_idx, -z_idx, +z_idx, mat_id] |


---

### Example

[geometry]

x -1.0 -0.5 0.0 0.5 1.0

y 0.0 1.0 63 lin

z 1e-3 1.0 16 log

fill 1

block 10 20 0 10 0 5 2