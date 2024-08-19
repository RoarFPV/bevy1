//! This example illustrates how to use [`States`] for high-level app control flow.
//! States are a powerful but intuitive tool for controlling which logic runs when.
//! You can have multiple independent states, and the [`OnEnter`] and [`OnExit`] schedules
//! can be used to great effect to ensure that you handle setup and teardown appropriately.
//!
//! In this case, we're transitioning from a `Menu` state to an `InGame` state.

use bevy::core_pipeline::tonemapping::Tonemapping;
use bevy::math::{UVec3, Vec3};
use bevy::utils::{HashSet, Hashed};
use bevy::{diagnostic::{FrameTimeDiagnosticsPlugin, EntityCountDiagnosticsPlugin}, pbr::CascadeShadowConfigBuilder, prelude::*};
use bevy_editor_pls::prelude::*;
use noisy_bevy::simplex_noise_2d_seeded;

use latticebolt::{
    CellDensity, CellDensityGenA, CellDensityGenB, CellId, LatticeConfig, LatticeRuntime,
};
// use rand::Rng;
use std::f32::consts::PI;

mod field;
mod latticebolt;
mod ui;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(FrameTimeDiagnosticsPlugin)
        .add_plugins(EntityCountDiagnosticsPlugin)
        //.add_plugins(EditorPlugin::new())
        .add_plugins(ui::UiPlugin)
        .register_type::<LatticeConfig>()
        //.register_type::<LatticeRuntime>()
        .register_type::<TerrainConfig>()
        .insert_resource(TerrainConfig {
            freq: 0.02,
            amp: 0.1,
            seed: 0,
        })
        .insert_resource(LatticeConfig {
            size: UVec3::new(100, 10, 50),
            rho0: 1.0,
            tau: 1.95,
            t_tau: 0.01,
            density_max: 0.4,
            view_density_max: 0.1,
            view_rho: 10.0,
            emitter: UVec3::new(30, 1, 30),
            ..default()
        })
        .add_systems(OnEnter(ui::State::Game), setup_game)
        .add_systems(
            Update,
            (run_game, camera_update, exit_app).run_if(in_state(ui::State::Game)),
        )
        // .add_systems(
        //     FixedUpdate,
        //     latticebolt::step.run_if(in_state(ui::State::Game)),
        // )
        .add_systems(
            FixedUpdate,
            (
                latticebolt::step_gen_a_to_b,
                latticebolt::step_gen_b_to_a,
                latticebolt::next_gen,
            )
                .chain()
                .run_if(in_state(ui::State::Game)),
        )
        .add_systems(OnExit(ui::State::Game), cleanup_game)
        .insert_resource(Time::<Fixed>::from_seconds(1.0 / 60.0))
        .run();
}

#[derive(Component, Debug)]
struct QuitButton;

#[derive(Resource)]
struct GameData {
    roots: HashSet<Entity>,
}

#[derive(Component, Debug)]
struct MainCamera;

#[derive(Default, Resource, Reflect)]
#[reflect(Resource)]
struct TerrainConfig {
    pub freq: f32,
    pub amp: f32,
    pub seed: u32,
}

fn setup_game(
    mut commands: Commands,
    config: Res<LatticeConfig>,
    terrain: Res<TerrainConfig>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    _asset_server: Res<AssetServer>,
) {
    let mut game = GameData {
        roots: HashSet::new(),
    };

    //commands.spawn(Camera2dBundle::default());

    let size = config.size;
    let camera = commands
        .spawn((
            Camera3dBundle {
                camera: Camera {
                    hdr: true,
                    ..default()
                },
                tonemapping: Tonemapping::AcesFitted,
                transform: Transform::from_xyz(
                    size.as_vec3().length(),
                    size.as_vec3().length(),
                    size.as_vec3().length() * 1.2,
                )
                .looking_at(size.as_vec3() / 2.0, Vec3::Y),
                ..default()
            },
            MainCamera,
        ))
        .id();

    game.roots.insert(
        commands
            .spawn(DirectionalLightBundle {
                directional_light: DirectionalLight {
                    shadows_enabled: false,
                    illuminance: 10000.0,
                    ..default()
                },
                transform: Transform {
                    translation: Vec3::new(0.0, 2.0, 0.0),
                    rotation: Quat::from_rotation_x(-PI / 4.),
                    ..default()
                },

                cascade_shadow_config: CascadeShadowConfigBuilder {
                    first_cascade_far_bound: 10.0,
                    maximum_distance: 200.0,
                    ..default()
                }
                .into(),
                ..default()
            })
            .id(),
    );

    game.roots.insert(camera);

    let mut rt = LatticeRuntime::default();

    let isize = config.size.as_ivec3();
    let center = isize / 2;

    let mut root = commands.spawn(SpatialBundle { ..default() });

    for z in 0..config.size.z {
        for x in 0..config.size.x {
            let n = simplex_noise_2d_seeded(
                Vec2::new(x as f32, z as f32) * terrain.freq,
                terrain.seed as f32,
            ) * terrain.amp;
            let h = (n.abs() * ((config.size.y * 4) as f32)) as u32;
            let mut hasbox = false;

            for y in 0..config.size.y {
                let pos = UVec3::new(x, y, z);
                let ipos = pos.as_ivec3();
                let mut solid = false;
                if h == y && !hasbox {
                    // game.roots.insert(commands.spawn(PbrBundle {
                    //     mesh: meshes.add(Mesh::from(shape::Box::new(1.0, 1.0, 1.0))),
                    //     material: materials.add(
                    //         Color::rgba_u8(
                    //             ((pos.y as f32 / config.size.y as f32) * 180.0) as u8,
                    //             230,
                    //             200,
                    //             255,
                    //         )
                    //         .into(),
                    //     ),
                    //     transform: Transform::from_translation(pos.as_vec3()),
                    //     ..default()
                    // }).id());
                    hasbox = true;
                    solid = true;
                }

                if y < h {
                    continue;
                }

                let index = LatticeRuntime::index(isize, ipos);
                let hash = Hashed::new(index);

                root.with_children(|parent| {
                    let child = parent
                        .spawn((
                            CellId {
                                id: index,
                                pos: ipos,
                            },
                            CellDensityGenA {
                                v: CellDensity {
                                    solid,
                                    density: latticebolt::density_init_eq(
                                        config.nodes,
                                        config.rho0,
                                        Vec3::X * config.density_max / 2.0,
                                    ),
                                    thermal_density: latticebolt::thermal_density_init_eq(
                                        config.nodes,
                                        config.rho0,
                                        Vec3::X * config.density_max / 2.0,
                                    ),
                                    ..default()
                                },
                                ..default()
                            },
                            CellDensityGenB {
                                v: CellDensity {
                                    solid,
                                    density: latticebolt::density_init_eq(
                                        config.nodes,
                                        config.rho0,
                                        Vec3::X * config.density_max / 2.0,
                                    ),
                                    thermal_density: latticebolt::thermal_density_init_eq(
                                        config.nodes,
                                        config.rho0,
                                        Vec3::X * config.density_max / 2.0,
                                    ),
                                    ..default()
                                },
                                ..default()
                            },
                        ))
                        .id();

                    rt.id_to_entity.insert(hash, child);
                });
            }
        }
    }

    game.roots.insert(root.id());
    // let id = latticebolt::LatticeRuntime::index(config.size.as_ivec3(), config.size.as_ivec3() / 2);
    // curr.cells[id].density[0] = config.density_max;

    commands.insert_resource(game);
    commands.insert_resource(rt);
}

fn cleanup_game(mut commands: Commands, mut game: ResMut<GameData>) {
    for entity in game.roots.iter() {
        commands.entity(*entity).despawn_recursive();
    }

    game.roots.clear();

    commands.remove_resource::<GameData>();
    commands.remove_resource::<LatticeRuntime>();
    // commands.remove_resource::<LatticeConfig>();
}

fn exit_app(keyboard: Res<ButtonInput<KeyCode>>, mut next_state: ResMut<NextState<ui::State>>) {
    if keyboard.just_pressed(KeyCode::Escape) {
        next_state.set(ui::State::Menu);
    }
}

fn run_game(
    lr: Res<LatticeRuntime>,
    lc: Res<LatticeConfig>,
    terrain: Res<TerrainConfig>,
    keyboard: Res<ButtonInput<KeyCode>>,
    mut cells: Query<(&CellId, &mut CellDensityGenA, &mut CellDensityGenB)>,

    mut gizmos: Gizmos,
) {
    for (id, a, b) in cells.iter() {
        let cell = if lr.gen { &b.v } else { &a.v };
        let normal_speed = (cell.total_velocity.length() / lc.density_max);

        //let density_scalar = cell.total_density.abs()/lc.view_density_max;

        // if density_scalar >= 0.01 {
        // gizmos.sphere(cell.pos.as_vec3(), Quat::IDENTITY, (1.0/(density_scalar*density_scalar))*0.5,
        //     Color::hsl(density_scalar * 360.0, 0.9, 0.5));
        // //((cell.total_density/lc.view_density_max + lc.view_density_max)/2.0);
        // }

        //print!("t: {:?}", cell.thermal_density[0]);
        gizmos.ray(
            id.pos.as_vec3(),
            (cell.total_velocity.normalize() * (normal_speed * normal_speed*2.0 + 0.1))
                .clamp_length_max(5.0), // * cell.density[0] / 5.0,
            Color::hsla(
                ((( cell.total_thermal_density/lc.view_density_max).clamp(-1.0, 1.0)+1.0)/2.0 * 360.0).clamp(180.0, 360.0),
                0.7,
                0.5 + if id.pos.x == 0 { 1.0 } else { 0.0 },
                0.1 + normal_speed / lc.view_rho, //(normal_speed).clamp(0.01, 1.0),
            ),
        );
    }

    if keyboard.pressed(KeyCode::Space) {
        let p = LatticeRuntime::index(lc.size.as_ivec3(), lc.emitter.as_ivec3());
        let hash = Hashed::new(p);

        if let Some(entity) = lr.id_to_entity.get(&hash) {
            if let Ok((_id, mut a, mut b)) = cells.get_mut(*entity) {
                let n = simplex_noise_2d_seeded(
                    lc.emitter.xy().as_vec2() * terrain.freq,
                    terrain.seed as f32,
                ) * terrain.amp;

                let change = lc.rho0 * n;

                a.v.density[0] = change; // [change; latticebolt::NODE_COUNT];
                b.v.density[0] = change; // = [change; latticebolt::NODE_COUNT];
            }
        }
    }

    // for cell in cells_next.cells.iter() {
    //     gizmos.ray(
    //         cell.pos.as_vec3(),
    //         Vec3::Z * (cell.density[0] / lc.rho0), //.clamp(-1.0, 1.0),
    //         Color::hsla((cell.density[0]) * 360.0 + 180.0, 0.7, 0.7, 0.5),
    //     );
    // }
}

fn camera_update(
    keyboard: Res<ButtonInput<KeyCode>>,
    mut q: Query<&mut Transform, With<MainCamera>>,
    lattice: Res<LatticeConfig>,
    time: Res<Time>,
) {
    let mut dir = Vec3::ZERO;

    if let Ok(mut transform) = q.get_single_mut() {
        if keyboard.pressed(KeyCode::KeyA) {
            dir.y += -1.0;
        }

        if keyboard.pressed(KeyCode::KeyD) {
            dir.y += 1.0;
        }

        if keyboard.pressed(KeyCode::ArrowUp) {
            dir.z += -1.0;
        }

        if keyboard.pressed(KeyCode::ArrowDown) {
            dir.z += 1.0;
        }

        if keyboard.pressed(KeyCode::KeyW) {
            dir.x += -1.0;
        }

        if keyboard.pressed(KeyCode::KeyS) {
            dir.x += 1.0;
        }

        if dir.length_squared() <= 0.0 {
            return;
        }

        dir *= time.delta_seconds();
        let forward = transform.local_z();
        let center = lattice.size.as_vec3() / 2.0;
        transform.translation += forward * dir.z * 10.0;
        transform.rotate_around(
            center,
            Quat::from_rotation_y(dir.y) * Quat::from_rotation_x(dir.x),
        );
        transform.look_at(center, Vec3::Y);
    }
}
