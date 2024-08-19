use bevy::prelude::*;

mod menu;

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq, Hash, States)]
pub enum State {
    #[default]
    Menu,
    Game,
}

#[derive(Resource)]
struct CameraData {
    camera: Entity,
}

const NORMAL_BUTTON: Color = Color::rgb(0.15, 0.15, 0.15);
const HOVERED_BUTTON: Color = Color::rgb(0.25, 0.25, 0.25);
const PRESSED_BUTTON: Color = Color::rgb(0.35, 0.75, 0.35);

pub struct UiPlugin;

impl Plugin for UiPlugin { 
    fn build(&self, app: &mut App){
        app
            .init_state::<State>()
            .add_plugins(menu::MenuPlugin)
            .add_systems(OnEnter(State::Menu), setup)
            .add_systems( OnExit(State::Menu), teardown)
            ;
    }
}

fn setup( mut commands: Commands) {
    let camera = commands.spawn(
        Camera2dBundle::default()
    ).id();

    commands.insert_resource(CameraData { camera });
}

fn teardown( mut commands: Commands, data: Res<CameraData>) {
   commands.entity(data.camera).despawn_recursive();
}