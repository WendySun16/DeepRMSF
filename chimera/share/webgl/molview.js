// Copyright © 2011 Regents of the University of California.
// All rights reserved.  This software provided pursuant to a
// license agreement containing restrictions on its disclosure,
// duplication and use.  This notice must be embedded in or
// attached to all copies, including partial copies, of the
// software or any revisions or derivations thereof.
//
'use strict';
// vi: sw=2:

function webGLStart() {
  var mouse_position;
  var canvas_name = 'molview';

  if (!PhiloGL.hasWebGL()) {
    return;
  }

  // from http://paulirish.com/2011/requestanimationframe-for-smart-animating/
  // shim layer with setTimeout fallback
  window.requestAnimFrame = (function() {
    return window.requestAnimationFrame
    || window.webkitRequestAnimationFrame
    || window.mozRequestAnimationFrame
    || window.oRequestAnimationFrame
    || window.msRequestAnimationFrame
    || function(/* function */ callback, /* DOMElement */ element) {
	 window.setTimeout(callback, 0);
       };
  })();

  // from http://www.jspatterns.com/category/patterns/code-reuse/
  function inherit(ChildClass, ParentClass) {
    var Chain = function () {};
    Chain.prototype = ParentClass.prototype;
    ChildClass.prototype = new Chain();
    ChildClass.prototype.constructor = ChildClass;
  }

  // grab some utility functions from PhiloGL
  function $() {};

  $.extend = function (to, from) {
    for (var p in from) {
      to[p] = from[p];
    }
    return to;
  };

  $.type = (function () {
    var oString = Object.prototype.toString,
	type = function (e) {
	  var t = oString.call(e);
	  return t.substr(8, t.length - 9).toLowerCase();
	};

    return function (elem) {
      var elemType = type(elem);
      if (elemType != 'object') {
	return elemType;
      }
      if (elem.$$family) return elem.$$family;
      return (elem && elem.nodeName && elem.nodeType == 1) ? 'element' : elemType;
    };
  })();

  (function () {
    function detach(elem) {
      var type = $.type(elem), ans;
      if (type == 'object') {
	ans = {};
	for (var p in elem) {
	  ans[p] = detach(elem[p]);
	}
	return ans;
      } else if (type == 'array') {
	ans = [];
	for (var i = 0, l = elem.length; i < l; i++) {
	  ans[i] = detach(elem[i]);
	}
	return ans;
      } else {
	return elem;
      }
    }

    $.merge = function () {
      var mix = {};
      for (var i = 0, l = arguments.length; i < l; i++){
	  var object = arguments[i];
	  if ($.type(object) != 'object') continue;
	  for (var key in object){
	      var op = object[key], mp = mix[key];
	      if (mp && $.type(op) == 'object' && $.type(mp) == 'object') {
		mix[key] = $.merge(mp, op);
	      } else{
		mix[key] = detach(op);
	      }
	  }
      }
      return mix;
    };
  })();

  PhiloGL.unpack();

  // create Spheres subclass of O3D.Model

  function Spheres(opt) {
    opt = $.merge({
      program: 'offset'
    }, opt || {});
    O3D.Model.call(this, opt);
    this.ProtoSpheres = {};
    this.spheres = {}
    this.render = Spheres.prototype.render;
  }

  inherit(Spheres, O3D.Model);

  $.extend(Spheres.prototype, {
      // disable normal object methods and do everything in render method
      setUniforms: function () {},
      setAttributes: function () {},
      setShininess: function () {},
      setReflection: function () {},
      setVertices: function () {},
      setColors: function () {},
      setPickingColors: function () {},
      setNormals: function () {},
      setTextures: function () {},
      setTexCoords: function () {},
      setIndices: function () {},
      unsetAttributes: function () {},
      unsetVertices: function () {},
      unsetColors: function () {},
      unsetPickingColors: function () {},
      unsetNormals: function () {},
      unsetTexCoords: function () {},
      unsetIndices: function () {},

      add: function (scene, radius, position, color) {
	var sphere = this.ProtoSpheres[radius];
	if (sphere == undefined) {
	  sphere = new O3D.Sphere({
	    radius: radius, nlat: 10, nlong: 10,
	    program: 'offset'
	  });
	  this.ProtoSpheres[radius] = sphere;
	  sphere.id = "sphere-" + radius;
	  scene.defineBuffers(sphere);
	}
	var data = this.spheres[radius];
	if (data == undefined) {
	  data = this.spheres[radius] = [];
	}
	data.push([position, color]);
      },

      render: function (gl, program, camera) {
	var offset = program.attributes['offset'];
	var color = program.attributes['color'];
	var obj;
	for (var radius in this.spheres) {
	  obj = this.ProtoSpheres[radius];
	  obj.setUniforms(program);
	  obj.setAttributes(program);
	  obj.setShininess(program);
	  obj.setReflection(program);
	  obj.setVertices(program);
	  obj.setColors(program);
	  obj.setPickingColors(program);
	  obj.setNormals(program);
	  obj.setTextures(program);
	  obj.setTexCoords(program);
	  obj.setIndices(program);

	  var data = this.spheres[radius];
	  for (var i = 0, l = data.length; i < l; ++i) {
	    gl.vertexAttrib3fv(offset, data[i][0]);
	    gl.vertexAttrib4fv(color, data[i][1]);
	    if (obj.indices) {
	      gl.drawElements((obj.drawType !== undefined)
			  ? gl.get(obj.drawType) : gl.TRIANGLES,
		  obj.indices.length, gl.UNSIGNED_SHORT, 0);
	    } else {
	      gl.drawArrays((obj.drawType !== undefined)
			  ? gl.get(obj.drawType) : gl.TRIANGLES,
		  0, obj.vertices.length / 3);
	    }
	  }

	  obj.unsetAttributes(program);
	  obj.unsetVertices(program);
	  obj.unsetColors(program);
	  obj.unsetPickingColors(program);
	  obj.unsetNormals(program);
	  obj.unsetTexCoords(program);
	  obj.unsetIndices(program);
	}
      },
  });

  // create Cylinder subclass of O3D.Model

  function Cylinders(opt) {
    opt = $.merge({
      program: 'cylinder'
    }, opt || {});
    O3D.Model.call(this, opt);
    this.ProtoCylinders = {};
    this.cylinders = {}
    this.render = Cylinders.prototype.render;
  }

  inherit(Cylinders, O3D.Model);

  $.extend(Cylinders.prototype, {
      // disable normal object methods and do everything in render method
      setUniforms: function () {},
      setAttributes: function () {},
      setShininess: function () {},
      setReflection: function () {},
      setVertices: function () {},
      setColors: function () {},
      setPickingColors: function () {},
      setNormals: function () {},
      setTextures: function () {},
      setTexCoords: function () {},
      setIndices: function () {},
      unsetAttributes: function () {},
      unsetVertices: function () {},
      unsetColors: function () {},
      unsetPickingColors: function () {},
      unsetNormals: function () {},
      unsetTexCoords: function () {},
      unsetIndices: function () {},

      add: function (scene, radius, height, mat4x3, color) {
	var cylinder = this.ProtoCylinders[radius];
	if (cylinder == undefined) {
	  cylinder = new O3D.Cylinder({
	    radius: radius, nvertical: 2, nradial: 10,
	    program: 'cylinder'
	  });
	  this.ProtoCylinders[radius] = cylinder;
	  cylinder.id = "cylinder-" + radius;
	  scene.defineBuffers(cylinder);
	}
	var data = this.cylinders[radius];
	if (data == undefined) {
	  data = this.cylinders[radius] = [];
	}
	data.push([height, mat4x3, color]);
      },

      render: function (gl, program, camera) {
	var yscale = program.attributes['yscale'];
	var transformX = program.attributes['transformX'];
	var transformY = program.attributes['transformY'];
	var transformZ = program.attributes['transformZ'];
	var color = program.attributes['color'];
	var obj;
	for (var radius in this.cylinders) {
	  obj = this.ProtoCylinders[radius];
	  obj.setUniforms(program);
	  obj.setAttributes(program);
	  obj.setShininess(program);
	  obj.setReflection(program);
	  obj.setVertices(program);
	  obj.setColors(program);
	  obj.setPickingColors(program);
	  obj.setNormals(program);
	  obj.setTextures(program);
	  obj.setTexCoords(program);
	  obj.setIndices(program);

	  var data = this.cylinders[radius];
	  for (var i = 0, l = data.length; i < l; ++i) {
	    gl.vertexAttrib1f(yscale, data[i][0]);
	    gl.vertexAttrib4fv(transformX, data[i][1][0]);
	    gl.vertexAttrib4fv(transformY, data[i][1][1]);
	    gl.vertexAttrib4fv(transformZ, data[i][1][2]);
	    gl.vertexAttrib4fv(color, data[i][2]);
	    if (obj.indices) {
	      gl.drawElements((obj.drawType !== undefined)
			  ? gl.get(obj.drawType) : gl.TRIANGLES,
		  obj.indices.length, gl.UNSIGNED_SHORT, 0);
	    } else {
	      gl.drawArrays((obj.drawType !== undefined)
			  ? gl.get(obj.drawType) : gl.TRIANGLES,
		  0, obj.vertices.length / 3);
	    }
	  }

	  obj.unsetAttributes(program);
	  obj.unsetVertices(program);
	  obj.unsetColors(program);
	  obj.unsetPickingColors(program);
	  obj.unsetNormals(program);
	  obj.unsetTexCoords(program);
	  obj.unsetIndices(program);
	}
      },
  });

  function buildScene(app, json) {
    var gl = app.gl;
    var camera = app.camera;
    var scene = app.scene;
    var lights = scene.config.lights;
    var seenLight = false;
    var index;
    var spheres, cylinders;
    var model;
	camera.wheelScale = 1;
    for (index in json) {
      var item = json[index];
      switch (item[0]) {
      case 's': // sphere
	if (!spheres) {
	  spheres = new Spheres;
	  scene.add(spheres);
	}
	spheres.add(scene, item[1], item[2], item[3]);
	break;
      case 'c': // cylinder
	if (!cylinders) {
	  cylinders = new Cylinders;
	  scene.add(cylinders);
	}
	cylinders.add(scene, item[1], item[2], item[3], item[4]);
	break;
      case 'p':
	model = new O3D.Model({
	  program: 'nolight',
	  drawType: "POINTS",
	  vertices: item[1],
	  colors: item[2]
	});
	scene.add(model);
	break;
      case 'l':
	model = new O3D.Model({
	  program: 'nolight',
	  drawType: "LINES",
	  vertices: item[1],
	  colors: item[2]
	});
	scene.add(model);
	break;
      case 'il':
	model = new O3D.Model({
	  program: 'nolight',
	  drawType: "LINES",
	  vertices: item[1],
	  colors: item[2],
	  indices: item[3]
	});
	scene.add(model);
	break;
      case 't':
      case 'ts':
	model = new O3D.Model({
	  program: 'default',
	  drawType: (item[0] == 't') ? "TRIANGLES" : "TRIANGLE_STRIP",
	  vertices: item[1],
	  normals: item[2],
	  colors: item[3],
	  indices: item[4]
	});
	scene.add(model);
	break;
      case 'bg': // background color
	gl.clearColor(item[1], item[2], item[3], 1);
	break;
      case 'vp': { // viewport
	  var canvas = app.canvas;
	  canvas.width = item[1];
	  canvas.height = item[2];
	  camera.near = item[3];
	  camera.far = item[4];
	  camera.aspect = item[1] / item[2];
	  camera.wheelScale = camera.near / 5;
	  if (camera.type == 'ortho') {
	    // figure out fov
	    camera.fov = Math.atan(camera.orthoParams[3] / camera.near) * 360 / Math.PI;
	  }
	};
	break;
      case 'la': // ambient light
	if (!seenLight) {
	  lights.enable = true;
	  lights.ambient = { r: 0, g: 0, b: 0 };
	  lights.directional = {
	    color: { r: 0, g: 0, b: 0 },
	    direction: { x: item[4], y: item[5], z: item[6] }
	  };
	  lights.points = [];
	  seenLight = true;
	};
	lights.ambient.r += item[1]
	lights.ambient.g += item[2]
	lights.ambient.b += item[3]
	break;
      case 'ld': // directional light
	if (!seenLight) {
	  lights.enable = true;
	  lights.ambient = { r: 0, g: 0, b: 0 };
	  lights.directional = {
	    color: { r: 0, g: 0, b: 0 },
	    direction: { x: item[4], y: item[5], z: item[6] }
	  };
	  lights.points = [];
	  seenLight = true;
	}
	var p = camera.target;
	var dir = new Vec3(-item[4], -item[5], -item[6]);
	dir.$scale(p.sub(camera.position).norm());
	lights.points.push({
	  diffuse: { r: item[1], g: item[2], b: item[3] },
	  specular: { r: 1.0, g: 1.0, b: 1.0 },
	  position: p.add(dir.scale(10))
	});
	break;
      case 'eyepos': // eye postion (look at eye position)
	camera.position = new Vec3(item[1], item[2], item[3]);
	break;
      case 'up': // up vector (look at up direction)
	camera.up = new Vec3(item[1], item[2], item[3]);
	break;
      case 'cofr': // center of rotation (look at point)
	camera.target = new Vec3(item[1], item[2], item[3]);
	break;
      case 'ortho': // orthographic viewpoint
	camera.type = 'ortho';
	camera.orthoParams = [item[1], item[2], item[3], item[4]]
	break;
      case 'persp': // perspective viewpoint
	camera.fov = item[1];
	break;
      }
    }
    camera.update();
  }

  function loadJSON(url, app) {
    new IO.XHR({
      url: url,
      onSuccess: function (text) {
	var json = JSON.parse(text);
	buildScene(app, json);
      }
    }).send();
    /*
    console.log('call JSONP');
    IO.JSONP({
      url: url,
      callbackKey: 'loadModel',
      onComplete: function(json) {
        console.log(json);
	buildScene(app, json);
      }
    });
    */
  }

  function draw(app) {
    var gl = app.gl,
	scene = app.scene,
	canvas = app.canvas;
    gl.viewport(0, 0, +canvas.width, +canvas.height);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    scene.render();
  }

  //Create application
  new PhiloGL(canvas_name, {
    program: [{
      id: 'default',
      from: 'defaults',
      vs: 'Default',
      fs: 'Default'
    },{
      id: 'offset',
      from: 'ids',
      vs: 'offset.vs',
      fs: 'default.fs',
      noCache: true		// TODO: false for production version
    },{
      id: 'cylinder',
      from: 'ids',
      vs: 'cylinder.vs',
      fs: 'default.fs',
    },{
      id: 'nolight',
      from: 'ids',
      vs: 'nolight.vs',
      fs: 'default.fs',
    }],
    /*
    context: {
      debug: true
    },
    */
    camera: {
      position: {
        x: 0, y: 0, z: -7
      }
    },
    events: {
      onDragStart: function (e) {
	var canvas = this.canvas;
	var radius = 0.9 * 0.5 * Math.min(canvas.width, canvas.height);
	mouse_position = {
	  x: e.x / radius,
	  y: e.y / radius
	};
      },
      onDragMove: function (e) {
	var canvas = this.canvas;
	var radius = 0.9 * 0.5 * Math.min(canvas.width, canvas.height);
	var new_position = {
	  x: e.x / radius,
	  y: e.y / radius
	};
	try {
	  var rotMat = vsphere(mouse_position, new_position);
	  var camera = this.camera;
	  var matrix = new Mat4;
	  matrix.$translate(camera.target[0], camera.target[1],
	    camera.target[2]);
	  matrix.$mulMat4(rotMat);
	  matrix.$translate(-camera.target[0], -camera.target[1],
	    -camera.target[2]);
	  for (var i = 0, models = this.scene.models, l = models.length; i < l; ++i) {
	    var elem = models[i];
	    elem.matrix = matrix.mulMat4(elem.matrix);
	  }
	} catch (err) {
	}
        mouse_position = new_position;
	var self = this;
	requestAnimFrame(function () { draw(self); });
      },
      onMouseWheel: function (e) {
        e.stop();
        var camera = this.camera;
	adjust = e.wheel * camera.wheelScale;
	if (camera.near + adjust > 0) {
	  camera.position.z += adjust;
	  camera.near += adjust;
	  camera.far += adjust;
	}
        camera.update();
	var self = this;
	requestAnimFrame(function () { draw(self); });
      }
    },
    onError: function (info) {
      alert("There was an error creating the app. " + info);
    },
    onLoad: function (app) {
      //Unpack app properties
      //app.gl = WebGLDebugUtils.makeDebugContext(app.gl);
      var gl = app.gl,
          scene = app.scene,
          canvas = app.canvas;

      var DEBUG_COUNT = 0;		// 0 to disable

      //Basic gl setup
      gl.clearColor(0.0, 0.0, 0.0, 1.0);
      gl.clearDepth(1.0);
      gl.enable(gl.DEPTH_TEST);
      gl.depthFunc(gl.LEQUAL);

      //Add objects to the scene
      //loadJSON('/gregc/one-cpk.json', app)
      //loadJSON('/gregc/mtx-cpk.json', app)
      //loadJSON('/gregc/3k9f-cpk.json', app)
      buildScene(app, json);

      requestAnimFrame(function () { draw(app); });
    }
  });
}
