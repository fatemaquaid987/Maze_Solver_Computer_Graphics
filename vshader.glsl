#version 120

attribute vec4 vPosition;
//attribute vec4 vColor;
//varying vec4 color;
attribute vec2 vTexCoord;
varying vec2 texCoord;

uniform mat4 ctm;
uniform mat4 projection;

void main()
{
	texCoord = vTexCoord;
	//color = vColor;
	gl_Position = projection * ctm * vPosition;
}
