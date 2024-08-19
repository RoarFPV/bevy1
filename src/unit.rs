use bevy::{
  math::{IVec3, UVec3, Vec3},
  prelude::*,
  utils::{petgraph::visit::IntoNeighbors, Hashed, PreHashMap},
};

pub struct UnitPlugin;

impl Plugin for UnitPlugin {
  fn build(&self, app: &mut App) {
      app.insert_resource(config);
  }
}


#[derive(Default, Resource, Reflect)]
#[reflect(Resource)]
struct UnitManager {
    
}



#[derive(Component, Debug)]
struct Unit {

}


#[derive(Component, Debug)]
struct Controller {

}
