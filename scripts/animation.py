# -*- coding: utf-8 -*-

from direct.showbase.ShowBase import ShowBase
from panda3d.core import WindowProperties
from panda3d.core import AmbientLight, DirectionalLight
from panda3d.core import Vec4

ShowBase().destroy()

class Game(ShowBase):
    
    def __init__(self):
        ShowBase.__init__(self)
        
        self.disableMouse()
        
        properties = WindowProperties()
        properties.setSize(1000, 750)
        self.win.requestProperties(properties)
        
        
        mainLight = DirectionalLight("main light")
        self.mainLightNodePath = render.attachNewNode(mainLight)
        self.mainLightNodePath.setHpr(45, -45, 0)
        render.setLight(self.mainLightNodePath)

        ambientLight = AmbientLight("ambient light")
        ambientLight.setColor(Vec4(0.2, 0.2, 0.2, 1))
        self.ambientLightNodePath = render.attachNewNode(ambientLight)
        render.setLight(self.ambientLightNodePath)
        
        render.setShaderAuto()

        self.camera.setPos(0, 0, 32)
        self.camera.setP(-90)



game = Game()
game.run()