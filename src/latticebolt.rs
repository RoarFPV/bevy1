use std::sync::Arc;

use bevy::{
    ecs::query::{BatchingStrategy, QueryEntityError},
    math::{IVec3, UVec3, Vec3},
    //utils::hashbrown::{hash_map, HashMap},
    prelude::*,
    //tasks::ComputeTaskPool,
    utils::{petgraph::visit::IntoNeighbors, Hashed, PreHashMap},
};

pub const NODE_COUNT: usize = 19;
pub const THERMAL_NODE_COUNT: usize = 7;

pub struct LaticeBoltPlugin;

impl Plugin for LaticeBoltPlugin {
    fn build(&self, app: &mut App) {
        let config = LatticeConfig::default();

        app.insert_resource(config);
    }
}

#[derive(Debug, Default, Copy, Clone, Reflect)]
pub struct Node {
    pub weight: f32,
    pub velocity: Vec3,
}

impl Node {
    fn new(weight: f32, velocity: Vec3) -> Self {
        Self { weight, velocity }
    }
}

#[derive(Resource, Debug, Reflect, Copy, Clone)]
#[reflect(Resource)]
pub struct LatticeConfig {
    // average density
    pub rho0: f32,
    // collision timescale
    pub tau: f32,
    pub t_tau: f32,
    pub density_max: f32,
    pub thermal_const: f32,
    pub view_rho: f32,
    pub view_density_max: f32,
    pub eq: f32,
    pub size: UVec3,
    pub emitter: UVec3,

    // lattice velocity and weights
    pub nodes: [Node; 19],
}

const REFLECTED_NODES: [usize; NODE_COUNT] = [
    0,  // 0 -> 0, center -> center
    3,  // 1 -> 3, north -> south
    4,  // 2 -> 4, east -> west
    1,  // 3 -> 1, south -> north
    2,  // 4 -> 2, west -> east
    6,  // 5 -> 6, top -> bottom
    5,  // 6 -> 5, bottom -> top
    9,  // 7 -> 9, north east -> south west
    10, // 8 -> 10, south east -> north west
    7,  // 9 -> 7, south west -> north east
    8,  // 10 -> 8, north west -> south east
    17, // 11 -> 17, top north -> bottom south
    18, // 12 -> 18, top east -> bottom west
    15, // 13 -> 15, top south -> bottom north
    16, // 14 -> 16, top west -> bottom east
    13, // 15 -> 13, bottom north -> top south
    14, // 16 -> 14, bottom east -> top west
    11, // 17 -> 11, bottom south -> top north
    12, // 18 -> 12, bottom west -> top east
];

impl Default for LatticeConfig {
    fn default() -> Self {
        Self {
            rho0: 100.0,
            view_rho: 0.0001,
            view_density_max: 100.0,
            tau: 1.8,
            t_tau: 1.1,
            thermal_const: 1.0,
            density_max: 1000.0,
            eq: 1.0,
            nodes: [
                Node::new(1.0 / 3.0, Vec3::new(0.0, 0.0, 0.0)),
                Node::new(1.0 / 18.0, Vec3::new(0.0, 1.0, 0.0)), // north
                Node::new(1.0 / 18.0, Vec3::new(1.0, 0.0, 0.0)), // east
                Node::new(1.0 / 18.0, Vec3::new(0.0, -1.0, 0.0)), // south
                Node::new(1.0 / 18.0, Vec3::new(-1.0, 0.0, 0.0)), // west
                Node::new(1.0 / 18.0, Vec3::new(0.0, 0.0, 1.0)), // top
                Node::new(1.0 / 18.0, Vec3::new(0.0, 0.0, -1.0)), // bottom
                Node::new(1.0 / 36.0, Vec3::new(1.0, 1.0, 0.0)), // north east
                Node::new(1.0 / 36.0, Vec3::new(1.0, -1.0, 0.0)), // south east
                Node::new(1.0 / 36.0, Vec3::new(-1.0, -1.0, 0.0)), // south west
                Node::new(1.0 / 36.0, Vec3::new(-1.0, 1.0, 0.0)), // north west
                Node::new(1.0 / 36.0, Vec3::new(0.0, 1.0, 1.0)), // top north
                Node::new(1.0 / 36.0, Vec3::new(1.0, 0.0, 1.0)), // top east
                Node::new(1.0 / 36.0, Vec3::new(0.0, -1.0, 1.0)), // top south
                Node::new(1.0 / 36.0, Vec3::new(-1.0, 0.0, 1.0)), // top west
                Node::new(1.0 / 36.0, Vec3::new(0.0, 1.0, -1.0)), // bottom north
                Node::new(1.0 / 36.0, Vec3::new(1.0, 0.0, -1.0)), // bottom east
                Node::new(1.0 / 36.0, Vec3::new(0.0, -1.0, -1.0)), // bottom south
                Node::new(1.0 / 36.0, Vec3::new(-1.0, 0.0, -1.0)), // bottom west
            ],
            size: UVec3::new(5, 5, 1),
            emitter: UVec3::ZERO,
        }
    }
}

#[derive(Default, Resource)]
pub struct LatticeRuntime {
    pub gen: bool,

    pub id_to_entity: PreHashMap<CellIndex, Entity>,
}

pub type CellIndex = u64;

impl LatticeRuntime {
    pub fn index(size: IVec3, pos: IVec3) -> CellIndex {
        // if pos.cmpge(size).any() || pos.cmplt(IVec3::ZERO).any() {
        //     return usize::MAX;
        // }

        if pos.x < 0
            || pos.x >= size.x
            || pos.y < 0
            || pos.y >= size.y
            || pos.z < 0
            || pos.z >= size.z
        {
            return CellIndex::MAX;
        }

        (pos.x + (size.x * pos.y) + (size.x * size.y * pos.z)) as CellIndex
    }

    pub fn _next_gen(&mut self) {
        self.gen = !self.gen;
    }
}

#[derive(Default, Copy, Clone)]
pub struct CellDensity {
    pub density: [f32; NODE_COUNT],
    pub total_density: f32,
    pub total_velocity: Vec3,
    pub solid: bool,

    pub total_thermal_density: f32,
    pub thermal_density: [f32; THERMAL_NODE_COUNT],
}

#[derive(Component, Default)]
// #[reflect(Component)]
pub struct CellDensityGenA {
    pub v: CellDensity,
}

#[derive(Component, Default)]
// #[reflect(Component)]
pub struct CellDensityGenB {
    pub v: CellDensity,
}

#[derive(Component, Default, Reflect, Copy, Clone)]
#[reflect(Component)]
pub struct CellId {
    pub id: CellIndex,
    pub pos: IVec3,
}

pub fn step_cell<F>(
    dt: f32,
    id_to_entity: &PreHashMap<CellIndex, Entity>,
    config: &LatticeConfig,
    read_id: &CellId,
    read: &CellDensity,
    write: &mut CellDensity,
    neighbors: F,
) where
    F: Fn(&Entity) -> Result<(CellId, CellDensity), QueryEntityError>,
{
    let mut next = *read;

    for i in 1..read.density.len() {
        let node = &config.nodes[i];
        // stream from outside
        let prev_index = LatticeRuntime::index(
            config.size.as_ivec3(),
            read_id.pos - node.velocity.as_ivec3(),
        );

        if let Some(neighbor_ent) = id_to_entity.get(&Hashed::new(prev_index)) {
            if let Ok((_id, n_density)) = neighbors(neighbor_ent) {
                next.density[i] = n_density.density[i];

                if i < THERMAL_NODE_COUNT {
                    next.thermal_density[i] = n_density.thermal_density[i];
                }
                //next.density[i] += (stream_cell.density[i] - next.density[i])*dt;
                // write.density[i] += (n_density.density[i]) * dt;

                // // make sure to account for density streamed out in this direction
                // write.density[i] -= (density) * dt;
                continue;
            }
        }
        // not found so reflect
        // let r = REFLECTED_NODES[i];
        next.density[i] = read.density[i];

        if i < THERMAL_NODE_COUNT {
            next.thermal_density[i] = read.thermal_density[i];
        }

        // let density = read.density[i];

        // let delta = (write.density[r]) * dt;
        // write.density[i] += delta;
        // write.density[i] -= density * dt;
    }

    let mut total_density: f32 = 0.0;
    let mut total_thermal_density: f32 = 0.0;
    for i in 0..next.density.len() {
        total_density += next.density[i];

        if i < THERMAL_NODE_COUNT {
            total_thermal_density += next.thermal_density[i];
        }
    }

    // update hydro dynamic velocity
    let mut vel = Vec3::ZERO;
    for (i, node) in config.nodes.iter().enumerate() {
        vel += node.velocity * next.density[i];
    }

    vel /= total_density;

    vel = vel.clamp_length_max(config.density_max);

    let p = read_id.pos;
    let mut w = dt;
    let mut tw = config.t_tau;
    if p.x == 0 {
        //vel = Vec3::new(config.density_max / 2.0, 0.0, 0.0);
        total_density = 1.0;
        total_thermal_density = 0.5;
        w = 1.0;
        tw = 1.0;
    } else if p.x == (config.size.x as i32 - 1) {
        //vel = Vec3::ZERO;
        total_density = 0.0;
        total_thermal_density = -0.5;
        w = 1.0;
        tw = 1.0;
    } else if p.y == 0 || p.z == 0{
        vel = Vec3::ZERO;
        total_density = 1.0;
        total_thermal_density = 0.0;
        w = 1.0;
        tw = 1.0;
    } else if p.y == config.size.y as i32 - 1 || p.z == config.size.z as i32 - 1{
        vel = Vec3::ZERO;
        total_density = 1.0;
        total_thermal_density = 0.0;
        w = 1.0;
        tw = 1.0;
    }

    if read.solid {
        vel = Vec3::ZERO;
        total_density = 1.0;
        total_thermal_density = 0.0;
        w = 1.0;
        tw=1.0;
    }


    if p == config.emitter.as_ivec3() {
        next.thermal_density[0] = config.eq;
        total_thermal_density = config.eq;
        total_density = 1.0;
        w = 1.0;
        tw = 1.0;
    }
    {
        // https://etheses.whiterose.ac.uk/13546/1/Thesis_for_web_new.pdf
        // LBM BGK model, where the temperature T (or the
        //     39
        //     internal energy e) is computed as a second order moment of the distribution functions f, as [62]

        //     and the temperature is related to the internal energy by e = (3/2) k/mT
        //     where k is the Boltzmann factor and m the mass [31].

        // Thermal lbm lattice
        // -  add second lattice D3Q7 for temperature
        // -  total_velocity used in thermal eqlib calc
        // -  thermal density[0] * 3e-4 added to +Y vel of density

        // Thermal equalization
        let v = vel;
        let velsq = vel.length_squared();

        for d in 0..THERMAL_NODE_COUNT {
            let n = config.nodes[d];

            let term = n.velocity.dot(v);
            let eqlib = n.weight
                * total_thermal_density
                * (1.0 + (3.0 * term) + (4.5) * term * term + (-1.5) * velsq);

            let density = next.thermal_density[d];

            let f = (1.0 - tw) * density + tw * eqlib;
            next.thermal_density[d] = f; //.max(0.0);
        }
    }

    next.total_thermal_density = total_thermal_density;
    // add any outside accelerations/velocities here
    //
    // Thermal velocity factor
    vel += Vec3::Y * next.thermal_density[0] * config.thermal_const;

    next.total_density = total_density;
    next.total_velocity = vel;

    let velsq = vel.length_squared();

    // collision function
    for d in 0..config.nodes.len() {
        let n = config.nodes[d];
        let term = n.velocity.dot(vel);
        let eqlib =
            n.weight * total_density * (1.0 + (3.0 * term) + (4.5) * term * term + (-1.5) * velsq);

        let density = next.density[d];

        let f = (1.0 - w) * density + w * eqlib;
        next.density[d] = f.max(0.0);
    }

    *write = next;
}

pub fn density_init_eq(nodes: [Node; NODE_COUNT], rho: f32, vel: Vec3) -> [f32; NODE_COUNT] {
    let mut density = [0.0; NODE_COUNT];

    let velsq = vel.length_squared();

    for d in 0..nodes.len() {
        let n = nodes[d];
        let term = n.velocity.dot(vel);
        density[d] = n.weight * rho * (1.0 + (3.0 * term) + (4.5) * term * term + (-1.5) * velsq);

        //write.density[d] = write.density[d].clamp(-config.density_max, config.density_max);
    }

    density
}

pub fn thermal_density_init_eq(
    nodes: [Node; NODE_COUNT],
    rho: f32,
    vel: Vec3,
) -> [f32; THERMAL_NODE_COUNT] {
    let mut density = [0.0; THERMAL_NODE_COUNT];

    let velsq = vel.length_squared();

    for d in 0..THERMAL_NODE_COUNT {
        let n = nodes[d];
        let term = n.velocity.dot(vel);
        density[d] = n.weight * rho * (1.0 + (3.0 * term) + (4.5) * term * term + (-1.5) * velsq);

        //write.density[d] = write.density[d].clamp(-config.density_max, config.density_max);
    }

    density
}

pub fn step_gen_a_to_b(
    config: Res<LatticeConfig>,
    runtime: ResMut<LatticeRuntime>,
    mut read: Query<(&CellId, &CellDensityGenA, &mut CellDensityGenB)>,
    neighbors: Query<(&CellId, &CellDensityGenA)>,
) {
    let gen = runtime.gen;
    if !gen {
        return;
    }

    let dt = config.tau;
    let nr = Arc::new(neighbors);
    read
        //.iter_mut()
        .par_iter_mut()
        .batching_strategy(BatchingStrategy::fixed(128))
        .for_each(|(id, read_gen, mut write_gen)| {
            let n = nr.clone();
            step_cell(
                dt,
                &runtime.id_to_entity,
                &config,
                id,
                &read_gen.v,
                &mut write_gen.v,
                move |e| {
                    if let Ok((id, cellgen)) = n.get(*e) {
                        Ok((*id, cellgen.v))
                    } else {
                        Err(QueryEntityError::NoSuchEntity(*e))
                    }
                },
            );
        });
}

pub fn step_gen_b_to_a(
    config: Res<LatticeConfig>,
    runtime: ResMut<LatticeRuntime>,
    mut read: Query<(&CellId, &CellDensityGenB, &mut CellDensityGenA)>,
    neighbors: Query<(&CellId, &CellDensityGenB)>,
) {
    let gen = runtime.gen;
    if gen {
        return;
    }

    let dt = config.tau;
    let nr = Arc::new(neighbors);
    read
        //.iter_mut()
        .par_iter_mut()
        .batching_strategy(BatchingStrategy::fixed(128))
        .for_each(|(id, read_gen, mut write_gen)| {
            let n = nr.clone();
            step_cell(
                dt,
                &runtime.id_to_entity,
                &config,
                id,
                &read_gen.v,
                &mut write_gen.v,
                move |e| {
                    if let Ok((id, cellgen)) = n.get(*e) {
                        Ok((*id, cellgen.v))
                    } else {
                        Err(QueryEntityError::NoSuchEntity(*e))
                    }
                },
            );
        });
}

pub fn next_gen(mut runtime: ResMut<LatticeRuntime>) {
    let gen = runtime.gen;
    runtime.gen = !gen;
}
