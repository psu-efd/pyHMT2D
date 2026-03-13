
import meshio
import vtk

import numpy as np

import shapefile

if __name__ == "__main__":

    reader = vtk.vtkPolyDataReader()
    reader.SetFileName("contour_line.vtk")
    reader.Update()

    polyData = reader.GetOutput()
    #print(polyData)

    points = np.array(polyData.GetPoints().GetData())
    print(points)

    lines = polyData.GetLines()
    print(lines)

    nLines = lines.GetNumberOfCells()
    print(nLines)

    line_data = polyData.GetLines().GetData()
    print(line_data)

    sourceIdList = vtk.vtkIdList()

    lines.InitTraversal()

    line_list = []

    #loop through all line segments
    cell_counter = 0
    while lines.GetNextCell(sourceIdList):
        print("cell_counter = ", cell_counter)

        pointCount = sourceIdList.GetNumberOfIds()

        idList = vtk.vtkIdList()

        current_list = []

        for idIndex in range(pointCount):
            sourceId = sourceIdList.GetId(idIndex)

            print(sourceId)

            current_list.append([points[sourceId][0], points[sourceId][1]])

        cell_counter += 1

        line_list.append(current_list)

    w = shapefile.Writer('line_test')
    w.field('name', 'C')

    w.line(line_list)

    w.record('linestring1')

    w.close()

    print("All done!")