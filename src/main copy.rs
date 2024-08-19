use bevy::prelude::*;
use bevy_editor_pls::prelude::*;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(EditorPlugin::default())
        .add_plugins(UnitPlugin)
        .run()
}

pub struct UnitPlugin;

impl Plugin for UnitPlugin { 
    fn build(&self, app: &mut App){
        app.add_systems(Startup, setup)
        .add_systems(Update, print_names)
        .add_systems(Update, print_actions);
    }
}

pub fn setup(mut commands: Commands) {
    commands.spawn((Unit { 
        name: "John".to_string(),
    }, Action { task: Task::Idle}));

    commands.spawn(Unit { 
        name: "2".to_string(),
    });
}

pub fn print_names( q: Query<&Unit, Without<Action>>) {
    println!("Units without actions:");
    for unit in q.iter() {
        println!("\tName: {}", unit.name);
    }
}

pub fn print_actions( q: Query<(&Unit, &Action), With<Action>>) {
    println!("Units with actions:");
    for (unit, action) in q.iter() {
        println!("Name: {}\n\taction: {}", unit.name, match action.task {
            Task::Idle => "Idle",
            Task::Move => "Move",
            Task::Attack => "Attack",
            Task::Die => "Die",
        });
    }
}

#[derive(Component)]
pub struct Unit {
    pub name: String,
}


#[derive(Component)]
pub struct Action { 
    pub task: Task,
}

#[derive(Debug)]
pub enum Task {
    Idle,
    Move,
    Attack,
    Die
}